#include "rtl_wfscann.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <fftw3.h>
#include <rtl-sdr.h>
#include <pthread.h>
#include <math.h>
#include <unistd.h>
#include <png.h>

//##################################main variables
struct globalsinfo{
	bool _continually = 0;
	bool _gnuplot = 0;
	bool _truncatefile = 0;
	int _output_format = _TXT_FORMAT_;
	char *_output_file = NULL;
	int _interval_seconds = -1;
}globals;
//##################################RTL variables
static rtlsdr_dev_t *rtl_dev;
struct rtlinfo{
	int _rate = 2048000;
	int _device = 0; 
	int _freq_start = 125000000; 
	int _freq_stop = 125000000; 
	int _gain = 0; 
	int _ppm = 0;
	uint8_t **_samplesBuffer;
	int *_currentFreq;
}rtl;
//##################################FFT variables
struct fftwInfo {
	int _bin_width  = 100;	
	int _window_type = _WINDOWS_Hanning;
	const double *_window;
	int _avg = 1;
	int _FFT_numSlices = 0;
	int _FFT_samplesPerSlice = 0;
	int _SPS_numSlices = 0;
	int _SPS_samplesPerSlice = 0;
	fftw_complex **_all_fftin;
	fftw_complex *_fftout;
	double **_pwr_fftw_out;
	double *_dB_fftw_out;
}fftw;
fftw_plan fftPlan = NULL;

//##################################threads
struct ThreadArgs {
    int start;
    int end;
    int FFT_samplesPerSlice;
    const uint8_t **samplesBuffer;
	fftw_complex **all_fftin;
};
pthread_t *threads = NULL;
struct ThreadArgs *threadArgsArray = NULL;
int numThreads = 0;
int iterationsPerThread = 0;
int remainingIterations = 0; 
//##################################RTL SDR
static int rtlsdr_init(){
	int device_count = rtlsdr_get_device_count();
	if (!device_count) {
		fprintf(stderr, "No supported devices found.\n");
		exit(1);
	}
	printf("#Starting rtl_map ~\n");
	printf("#Found %d device(s):\n", device_count);
	for(int i = 0; i < device_count; i++){
		printf("#%d: %s\n", i, rtlsdr_get_device_name(i));
	}
	int dev_open = rtlsdr_open(&rtl_dev, rtl._device);
	if (dev_open < 0) {
		fprintf(stderr, "Failed to open RTL-SDR device #%d\n", rtl._device);
		exit(1);
	}else{
		printf("#Using device: #%d\n", dev_open);
	}
	if(!rtl._gain){
		rtlsdr_set_tuner_gain_mode(rtl_dev, rtl._gain);
		printf("#Gain mode set to auto.\n");
	}else{
		rtlsdr_set_tuner_gain_mode(rtl_dev, 1);
		int gain_count = rtlsdr_get_tuner_gains(rtl_dev, NULL);
		printf("#Supported gain values (%d): ", gain_count);
		int gains[gain_count], supported_gains = rtlsdr_get_tuner_gains(rtl_dev, gains), target_gain=0;
		for (int i = 0; i < supported_gains; i++){
			if (gains[i] < rtl._gain*10){target_gain = gains[i];}
		}
		printf("\n");
		printf("#Gain set to %.1f\n", target_gain / 10.0);
		rtlsdr_set_tuner_gain(rtl_dev, target_gain);
	}
	rtlsdr_set_freq_correction(rtl_dev, rtl._ppm);
	rtlsdr_set_center_freq(rtl_dev, rtl._freq_start);
	rtlsdr_set_sample_rate(rtl_dev, rtl._rate);
	printf("#Frequency start set to %d Hz.\n", rtl._freq_start);
	printf("#Sampling at %d sps\n", rtl._rate);
	int r = rtlsdr_reset_buffer(rtl_dev);
	if (r < 0){
		fprintf(stderr, "Failed to reset buffers.\n");
		return 1;
	}
	return 0;
}

//##################################function for images
struct ImageData{
    int _width;
	int _height = 1;
	int _actualLine;
    uint8_t** _data; // Tableau 1D de pixels (4 canaux par pixel)
	uint8_t _compress;
	int _max_value_db;
	int _min_value_db;
	double _scaleFttToPng; //255.0 / (max_value_db - min_value_db)
} Image;

void Init_image_data(int width, int height = 1, uint8_t compress = 0, int max_value_db = 0, int min_value_db = -110)
{
	Image._width = width;
	Image._height = height;
	Image._compress = compress;
	Image._max_value_db = max_value_db;
	Image._min_value_db = min_value_db;
	Image._actualLine = 0;
	Image._scaleFttToPng = 255.0 / (max_value_db - (min_value_db));
	Image._data = (uint8_t**)malloc(height * sizeof(uint8_t *));
	for (int i = 0; i < height; ++i){
		Image._data[i] = (uint8_t*)calloc(width * 4, sizeof(uint8_t));
	}
}

void convert_db_data_to_image_png(FILE *fp,double* fftw_data) {
	
	//prep data
	if(Image._data[Image._actualLine]){
		Image._data[Image._actualLine] = (uint8_t*)malloc(Image._width * 4);
	}
    for (int i = 0; i < Image._width; ++i) {
        int index = (int)(((fftw_data[i] - Image._min_value_db) * Image._scaleFttToPng));
        index = (index < 0) ? 0 : ((index > 255) ? 255 : index);
        int base_index = i * 4;
        Image._data[Image._actualLine][base_index] = websdrWaterfall[index][0];
        Image._data[Image._actualLine][base_index + 1] = websdrWaterfall[index][1];
        Image._data[Image._actualLine][base_index + 2] = websdrWaterfall[index][2];
        Image._data[Image._actualLine][base_index + 3] = 255; // Alpha
    }
	Image._actualLine = (Image._actualLine < Image._height) ? (Image._actualLine + 1) : 0;	
	
	//write data
	if (!fp) {
        fprintf(stderr, "Error opening file for writing.\n");
        return;
    }
    png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png) {
        fprintf(stderr, "Error creating PNG structure.\n");
        return;
    }
    png_infop info = png_create_info_struct(png);
    if (!info) {
        fprintf(stderr, "Error creating PNG information structure.\n");
        png_destroy_write_struct(&png, NULL);
        return;
    }
    if (setjmp(png_jmpbuf(png))) {
        fprintf(stderr, "Error initializing data for writing PNG.\n");
        png_destroy_write_struct(&png, &info);
        return;
    }
    png_init_io(png, fp);
    png_set_IHDR(
        png,
        info,
        Image._width, Image._height,
        8, PNG_COLOR_TYPE_RGBA, PNG_INTERLACE_NONE,
        PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT
    );
    png_write_info(png, info);
	int index=0;
    for (int y = Image._actualLine; y < Image._height; y++) {
        png_write_row(png, Image._data[index]);
		index++;
    }
    for (int y = 0; y < Image._actualLine-1; y++) {
        png_write_row(png, Image._data[index]);
		index++;
    }
	
    png_write_end(png, NULL);
    png_destroy_write_struct(&png, &info);
}


//##################################function for options
void help() {
    printf("Usage: rtl_wfscann [OPTIONS]\n");
    printf("Options:\n");
    printf("  -h, --help               Display this help message\n");
    printf("  -c, --continuous         Run continuously\n");
    printf("  -d, --device DEVICE_ID   Set RTL device ID\n");
    printf("  -f, --freq FREQ_RANGE    Set frequency range option (e.g., 125000000:155000000:100)\n");
    printf("  -g, --gain GAIN          Set RTL gain, max 50\n");
    printf("  -p, --ppm PPM            Set RTL PPM\n");
    printf("  -r, --rate RATE          Set RTL rate\n");
    printf("  -w, --window WINDOW      Set FFT window\n");
    printf("  -a, --avg AVG            Set FFT averaging\n");
    printf("  -G, --gnuplot            Include Gnuplot option and force FORMAT to TXT (e.g., plot results)\n");
    printf("  -O, --output FILE        Set output file\n");
    printf("  -I, --interval SECONDS   Set measurement interval in seconds 0 or -1 is disable as default\n");
    printf("  -F, --format FORMAT      Set output format (TXT, IMG, RAW) default is TXT\n");
    printf("  -T, --truncatefile       Set truncated output filename like timestamp_filename.txt\n");    
	printf("\nExample:\n");
    printf("./rtl_wfscann -f 143000000:147000000:2000 -g 10\n");
    printf("./rtl_wfscann -f 143000000:147000000:2000 -g 10 -c -G -o output.txt -i 5 -F TXT\n");
}

void validateFreqOption(char *freqRange) {
    char *token = strtok(freqRange, ":");

    if (token != NULL) {
        rtl._freq_start = atoi(token);
    } else {
        fprintf(stderr, "#Error: Missing value for --freq\n");
        help();
        exit(1);
    }

    token = strtok(NULL, ":");

    if (token != NULL) {
        rtl._freq_stop = atoi(token);
    } else {
        fprintf(stderr, "#Error: Missing value for --freq\n");
        help();
        exit(1);
    }

    token = strtok(NULL, ":");

    if (token != NULL) {
        fftw._bin_width = atoi(token);
    } else {
        fprintf(stderr, "#Error: Missing value for --freq\n");
        help();
        exit(1);
    }

    if (rtl._freq_start < 0 || rtl._freq_stop < 0 || fftw._bin_width <= 0 || rtl._freq_start >= rtl._freq_stop) {
        fprintf(stderr, "#Error: Values for --freq are invalid\n");
        help();
        exit(1);
    }

	int adjustedBandwidth = 0;
	while((rtl._freq_start+adjustedBandwidth*rtl._rate)<rtl._freq_stop)adjustedBandwidth++;
	rtl._freq_stop = rtl._freq_start+adjustedBandwidth*rtl._rate;

    printf("#Frequency start : %d\n", rtl._freq_start);
    printf("#Frequency stop  : %d\n", rtl._freq_stop);
    printf("#FFT binary width: %d\n", fftw._bin_width);
}

int validateFormatOption(const char *format) {
    if (strcmp(format, "TXT") == 0) {
        return _TXT_FORMAT_;
    } else if (strcmp(format, "IMG") == 0) {
        return _IMG_FORMAT_;
    } else if (strcmp(format, "RAW") == 0) {
        return _RAW_FORMAT_;
    } else {
        fprintf(stderr, "Unsupported output format: %s (TXT|IMG|RAW)\n", format);
        help();
        exit(1);
    }
}

int validateFFTWindowingOption(const char *window) {
    if (strcmp(window, "rectangle") == 0) {
        return _WINDOWS_Rectangle;
    } else if (strcmp(window, "hanning") == 0) {
        return _WINDOWS_Hanning;
    } else if (strcmp(window, "hamming") == 0) {
        return _WINDOWS_Hamming;
    } else if (strcmp(window, "blackman") == 0) {
        return _WINDOWS_Blackman;
    } else if (strcmp(window, "triangulaire") == 0) {
        return _WINDOWS_Triangulaire;
    } else if (strcmp(window, "bartlett") == 0) {
        return _WINDOWS_Bartlett;
    } else {
        fprintf(stderr, "Unsupported FFT windowing option: %s (rectangle|hanning|hamming|blackman|triangulaire|bartlett)\n", window);
        exit(1);
    }
}

char* generateFileNameWithTimestamp(const char* _output_file) {
    if (_output_file == NULL) {
        return NULL;
    }
    time_t timestamp = time(NULL);
    struct tm* timeinfo = localtime(&timestamp);
    if (timeinfo == NULL) {
        return NULL;
    }
    char timestampStr[20];
    strftime(timestampStr, sizeof(timestampStr), "%Y%m%d%H%M%S", timeinfo);
    size_t newFileNameLen = strlen(_output_file) + 20 + 2;
    char* newFileName = (char*)malloc(newFileNameLen);
    if (newFileName == NULL) {
        return NULL;
    }
    const char* lastSlash = strrchr(_output_file, '/');
    if (lastSlash != NULL) {
        size_t prefixLen = lastSlash - _output_file + 1;
        strncpy(newFileName, _output_file, prefixLen);
        sprintf(newFileName + prefixLen, "%s_%s", timestampStr, lastSlash + 1);
    } else {
        sprintf(newFileName, "%s_%s", timestampStr, _output_file);
    }
    return newFileName;
}
//##################################main
const double *generateFFTWindow(int length, int windowType) {
    double *window = (double *)malloc(length * sizeof(double));
    if (window == NULL) {
        fprintf(stderr, "Memory allocation error for FFT window.\n");
        exit(1);
    }

    switch (windowType) {
        case _WINDOWS_Rectangle:
            for (int i = 0; i < length; i++) {
                window[i] = 1.0;
            }
            break;

        case _WINDOWS_Hanning:
            for (int i = 0; i < length; i++) {
                window[i] = 0.5 * (1.0 - cos(2.0 * M_PI * i / (length - 1)));
            }
            break;

        case _WINDOWS_Hamming:
            for (int i = 0; i < length; i++) {
                window[i] = 0.54 - 0.46 * cos(2.0 * M_PI * i / (length - 1));
            }
            break;

        case _WINDOWS_Blackman:
            for (int i = 0; i < length; i++) {
                window[i] = 0.42 - 0.5 * cos(2.0 * M_PI * i / (length - 1)) + 0.08 * cos(4.0 * M_PI * i / (length - 1));
            }
            break;

        case _WINDOWS_Triangulaire:
            for (int i = 0; i < length; i++) {
                window[i] = 1.0 - fabs((2.0 * i / (length - 1)) - 1.0);
            }
            break;

        case _WINDOWS_Bartlett:
            for (int i = 0; i < length; i++) {
                window[i] = 1.0 - fabs((2.0 * i / (length - 1)) - 1.0);
            }
            break;

        default:
            fprintf(stderr, "Unsupported window type.\n");
            exit(1);
    }

    return window;
}

void* processPrepSamples(void* args) {
    struct ThreadArgs* threadArgs = (struct ThreadArgs*)args;
    int FFT_samplesPerSlice = threadArgs->FFT_samplesPerSlice;
    for (int i = threadArgs->start; i < threadArgs->end; i++) {
        for (int j = 0; j < FFT_samplesPerSlice; j++) {
			threadArgs->all_fftin[i][j][_Q_] = (threadArgs->samplesBuffer[i][j*2] - 127.4) * (1.0 / 128.0) * fftw._window[j];
			threadArgs->all_fftin[i][j][_I_] = (threadArgs->samplesBuffer[i][j*2+1] - 127.4) * (1.0 / 128.0) * fftw._window[j];
        }
    }

    pthread_exit(NULL);
}

#include <sys/time.h>
struct timeval start_time, end_time;
void ReadAndRun(){

	// //////////////////
	// gettimeofday(&start_time, NULL);
	// //////////////////
		
	int n_read=0;
	for (int j = 0; j < fftw._SPS_numSlices; j=j+fftw._avg) {
		rtlsdr_set_center_freq(rtl_dev, rtl._currentFreq[j]);
		rtlsdr_reset_buffer(rtl_dev);
		for (int i = 0; i < fftw._avg; i++) {
			if (rtlsdr_read_sync(rtl_dev, rtl._samplesBuffer[j+i], fftw._SPS_samplesPerSlice, &n_read) < 0) {
				fprintf(stderr, "Memory allocation error for sample buffer.\n");		
			}
			if(n_read!=fftw._SPS_samplesPerSlice){
				fprintf(stderr, "Error reading samples %d/%d.\n",n_read,fftw._SPS_samplesPerSlice);		
			}
		}
	}

	// ////////////////////////////////////
	// gettimeofday(&end_time, NULL);
	// long elapsed_time = (end_time.tv_sec - start_time.tv_sec) * 1000000L + (end_time.tv_usec - start_time.tv_usec);
	// double elapsed_time_seconds = (double)elapsed_time / 1000000.0;
	// printf("#Execution time : %.6f microsecondes\n", elapsed_time_seconds);
	// ////////////////////////////////////
	// //////////////////
	// gettimeofday(&start_time, NULL);
	// //////////////////
	
	iterationsPerThread = fftw._SPS_numSlices / numThreads;
	remainingIterations = fftw._SPS_numSlices % numThreads;  
	
	int start = 0;
    for (int i = 0; i < numThreads; i++) {
        struct ThreadArgs *args = &threadArgsArray[i];
        args->start = start;
        args->end = start + iterationsPerThread + (remainingIterations > 0 ? 1 : 0);
		remainingIterations--;
		
		args->FFT_samplesPerSlice = fftw._FFT_samplesPerSlice;
		args->samplesBuffer = (const uint8_t **)rtl._samplesBuffer;
		args->all_fftin = fftw._all_fftin;

        pthread_create(&threads[i], NULL, processPrepSamples, (void*)args);
		start = args->end;
    }

	for (int i = 0; i < fftw._FFT_numSlices; i++) {
			memset(fftw._pwr_fftw_out[i], 0, fftw._FFT_samplesPerSlice * sizeof(double));
	}	
		
    for (int i = 0; i < numThreads; i++) {
        pthread_join(threads[i], NULL);
    }

	// ////////////////////////////////////
	// gettimeofday(&end_time, NULL);
	// elapsed_time = (end_time.tv_sec - start_time.tv_sec) * 1000000L + (end_time.tv_usec - start_time.tv_usec);
	// elapsed_time_seconds = (double)elapsed_time / 1000000.0;
	// printf("#Execution time : %.6f \n", elapsed_time_seconds);
	// ////////////////////////////////////
	// //////////////////
	// gettimeofday(&start_time, NULL);
	// //////////////////
	
	for (int i = 0; i < fftw._SPS_numSlices; i++) {
		fftw_execute_dft(fftPlan, fftw._all_fftin[i], fftw._fftout);
		int index = i / fftw._avg;
		int k = 0;
		for (int j = 0; j < fftw._FFT_samplesPerSlice; j++) {
			if(j<(fftw._FFT_samplesPerSlice/2)) {k=(fftw._FFT_samplesPerSlice/2)-j;}
			else{k=(fftw._FFT_samplesPerSlice)-(j-(fftw._FFT_samplesPerSlice/2));}
			fftw._pwr_fftw_out[index][j] += fftw._fftout[k][_Q_] * fftw._fftout[k][_Q_] + fftw._fftout[k][_I_] * fftw._fftout[k][_I_];
		}
	}
	
	// ////////////////////////////////////
	// gettimeofday(&end_time, NULL);
	// elapsed_time = (end_time.tv_sec - start_time.tv_sec) * 1000000L + (end_time.tv_usec - start_time.tv_usec);
	// elapsed_time_seconds = (double)elapsed_time / 1000000.0;
	// printf("#Execution time : %.6f \n", elapsed_time_seconds);
	// ////////////////////////////////////
	// //////////////////
	// gettimeofday(&start_time, NULL);
	// //////////////////
	
	for (int i = 0; i < fftw._FFT_numSlices; i++) {
		for (int j = 0; j < fftw._FFT_samplesPerSlice; j++) {
			fftw._dB_fftw_out[i*fftw._FFT_samplesPerSlice + j] = round(((10*log10(fftw._pwr_fftw_out[i][j]/fftw._avg))+ db_adc_rtl_factor)* 100) / 100;
		}
	}
	
	// ////////////////////////////////////
	// gettimeofday(&end_time, NULL);
	// elapsed_time = (end_time.tv_sec - start_time.tv_sec) * 1000000L + (end_time.tv_usec - start_time.tv_usec);
	// elapsed_time_seconds = (double)elapsed_time / 1000000.0;
	// printf("#Execution time : %.6f \n", elapsed_time_seconds);
	// ////////////////////////////////////
	// //////////////////
	// gettimeofday(&start_time, NULL);
	// //////////////////
	
	if(globals._output_file) {
		FILE *fichier = NULL;
		if(globals._truncatefile){
			fichier = fopen(generateFileNameWithTimestamp(globals._output_file), "w");
		}else{
			fichier = fopen(globals._output_file, "w");
			}
		switch (globals._output_format) {
			case _TXT_FORMAT_:
				for (int i = 0; i < fftw._FFT_numSlices*fftw._FFT_samplesPerSlice; i++) {
					fprintf(fichier, "%d %f\n", rtl._freq_start + (i * fftw._bin_width), fftw._dB_fftw_out[i]);
				}
				break;
			case _IMG_FORMAT_:
				convert_db_data_to_image_png(fichier,fftw._dB_fftw_out);
				break;
			case _RAW_FORMAT_:
				fwrite(fftw._dB_fftw_out, sizeof(double), fftw._FFT_numSlices * fftw._FFT_samplesPerSlice, fichier);
				break;
			default:
				break;
			}
		fclose(fichier);
	}
	else{
		switch (globals._output_format) {
			case _TXT_FORMAT_:
				if(globals._gnuplot){printf("e\nplot \"-\"\n");}
				for (int i = 0; i < fftw._FFT_numSlices*fftw._FFT_samplesPerSlice; i++) {
					printf("%d %f\n", rtl._freq_start + (i * fftw._bin_width), fftw._dB_fftw_out[i]); 
				}
				break;
			case _RAW_FORMAT_:
				fwrite(fftw._dB_fftw_out, sizeof(double), fftw._FFT_numSlices * fftw._FFT_samplesPerSlice, stdout);
				break;
			default:
				break;
			}
	}
	// ////////////////////////////////////
	// gettimeofday(&end_time, NULL);
	// elapsed_time = (end_time.tv_sec - start_time.tv_sec) * 1000000L + (end_time.tv_usec - start_time.tv_usec);
	// elapsed_time_seconds = (double)elapsed_time / 1000000.0;
	// printf("#Execution time : %.6f \n", elapsed_time_seconds);
	// ////////////////////////////////////
	
}

void VarAndBuffInit(){
	
	fftPlan = fftw_plan_dft_1d((rtl._rate / fftw._bin_width), NULL, NULL, FFTW_BACKWARD, FFTW_ESTIMATE);

	fftw._FFT_numSlices = (rtl._freq_stop - rtl._freq_start) / rtl._rate;
	fftw._FFT_samplesPerSlice = (rtl._rate / fftw._bin_width); // /2 because fftw_complex is I/Q
	fftw._SPS_numSlices = fftw._FFT_numSlices * fftw._avg;
	fftw._SPS_samplesPerSlice = fftw._FFT_samplesPerSlice * 2;// *2 beacause 1 sample = 2 byte I/Q
	
	rtl._currentFreq = (int *)malloc(fftw._SPS_numSlices * sizeof(int));
	if (rtl._currentFreq == NULL) {
		fprintf(stderr, "Memory allocation error for frequencies index.\n");
		exit(1);
	}
	for (int j = 0; j < fftw._SPS_numSlices; j=j+fftw._avg) {
		rtl._currentFreq[j] = rtl._freq_start + (j/fftw._avg * rtl._rate) + (rtl._rate / 2);
	}
	
	fftw._window = generateFFTWindow(fftw._FFT_samplesPerSlice,fftw._window_type);	

	rtl._samplesBuffer = (uint8_t **)malloc(fftw._SPS_numSlices * sizeof(uint8_t *));
    if (rtl._samplesBuffer == NULL) {
        fprintf(stderr, "Memory allocation error for sample buffer.\n");
        exit(1);
    }
	
	for (int i = 0; i < fftw._SPS_numSlices; i++) {
		rtl._samplesBuffer[i] = (uint8_t *)malloc(fftw._SPS_samplesPerSlice * sizeof(uint8_t));
        if (rtl._samplesBuffer[i] == NULL) {
            fprintf(stderr, "Memory allocation error for sample buffer.\n");
            exit(1);
        }
	}
	
	fftw._all_fftin = (fftw_complex **)fftw_malloc(fftw._SPS_numSlices * sizeof(fftw_complex *));
	if (fftw._all_fftin == NULL) {
        fprintf(stderr, "Memory allocation error for fft input buffer.\n");
        exit(1);
    }
	for (int i = 0; i < fftw._SPS_numSlices; i++) {
		fftw._all_fftin[i]= (fftw_complex *)fftw_malloc(fftw._FFT_samplesPerSlice * sizeof(fftw_complex));
        if (fftw._all_fftin[i] == NULL) {
            fprintf(stderr, "Memory allocation error for fft input buffer.\n");
            exit(1);
        }
	}
	
	fftw._fftout = (fftw_complex *)fftw_malloc(fftw._FFT_samplesPerSlice * sizeof(fftw_complex));
    if (fftw._fftout == NULL) {
        fprintf(stderr, "Memory allocation error for fft output buffer\n");
        exit(1);
    }
	
	fftw._pwr_fftw_out = (double **)malloc(fftw._FFT_numSlices * sizeof(double *));
    if (fftw._pwr_fftw_out == NULL) {
        fprintf(stderr, "Memory allocation error for fft power values.\n");
        exit(1);
    }
	for (int i = 0; i < fftw._FFT_numSlices; i++) {
		fftw._pwr_fftw_out[i] = (double *)calloc(fftw._FFT_samplesPerSlice,sizeof(double));
		if (fftw._pwr_fftw_out[i] == NULL) {
			fprintf(stderr, "Memory allocation error for fft power values.\n");
			exit(1);
		}
	}
	
	numThreads = sysconf(_SC_NPROCESSORS_ONLN);
	threads = (pthread_t*)malloc(numThreads * sizeof(pthread_t));
	threadArgsArray = (struct ThreadArgs*)malloc(numThreads * sizeof(struct ThreadArgs));
	if (threadArgsArray == NULL) {
		fprintf(stderr, "Thread allocation error for prep threads.\n");
		exit(1);
	}
	
	fftw._dB_fftw_out = (double *)calloc(fftw._FFT_numSlices*fftw._FFT_samplesPerSlice,sizeof(double));
	if (fftw._dB_fftw_out == NULL) {
		fprintf(stderr, "Memory allocation error for fft power values in dB.\n");
		exit(1);
	}
	
	if(globals._output_file) {
		
		char* new_filename = (char*)malloc(strlen(globals._output_file) + strlen(".info") + 1);
		strcpy(new_filename, globals._output_file);strcat(new_filename, ".info");
		
		FILE *fichier = fopen(new_filename, "w");
		printf("#File name: %s\n", globals._output_file);
		if (fichier == NULL) {
			fprintf(stderr, "Impossible to open the file.\n");
			exit(1);
		}
		fprintf(fichier, "freq_start:%d\nfreq_stop:%d\nbin_bandwith:%d\ntotal_samples:%d\n", rtl._freq_start, rtl._freq_stop,fftw._bin_width,fftw._FFT_numSlices*fftw._FFT_samplesPerSlice);
		fclose(fichier);
		
		if(globals._output_format == _IMG_FORMAT_){
			if(globals._truncatefile){
			//  Init_image_data(int width , int height = 1, uint8_t compress = 0, int max_value_db = 0, int min_value_db = -110)
				Init_image_data(fftw._FFT_numSlices*fftw._FFT_samplesPerSlice, 1    , 0, 0, -110);
			}else{
				Init_image_data(fftw._FFT_numSlices*fftw._FFT_samplesPerSlice, 86400, 0, 0, -110);
			}
		}
	}

}

void VarAndBuffDeInit(){
	
	free(threadArgsArray);
	free(threads);
	
	for (int i = 0; i < fftw._FFT_numSlices; i++) {
		free(fftw._pwr_fftw_out[i]);
	}
	free(fftw._pwr_fftw_out);
	
	fftw_free(fftw._fftout);
	
	for (int i = 0; i < fftw._SPS_numSlices; i++) {
		fftw_free(fftw._all_fftin[i]);
	}
	fftw_free(fftw._all_fftin);

	for (int i = 0; i < fftw._SPS_numSlices; i++) {
		free(rtl._samplesBuffer[i]);
	}
	free(rtl._samplesBuffer);
	
	free((void *)fftw._window);
	fftw_destroy_plan(fftPlan);
	fftw_cleanup();
}

int main(int argc, char * argv[]) {
    int x;
	setbuf(stdout, NULL);
	
	static struct option long_options[] =
        {
			{"help",  no_argument, 0, 'h'},
			{"continuous",  no_argument,       0, 'c'},
			{"device",  required_argument, 0, 'd'},
			{"freq",    required_argument, 0, 'f'},
			{"gain",    required_argument, 0, 'g'},
			{"ppm",    required_argument, 0, 'p'},
			{"rate",    required_argument, 0, 'r'},
			{"window",    required_argument, 0, 'w'},
			{"avg",    required_argument, 0, 'a'},
			{"gnuplot", no_argument, 0, 'G'},
			{"output", required_argument, 0, 'O'},
			{"interval", required_argument, 0, 'I'},
			{"format", required_argument, 0, 'F'},
			{"truncatefile", no_argument, 0, 'T'},
			{0, 0, 0, 0}
        };
	int option_index = 0;

    while ((x = getopt_long(argc, argv, "hcd:f:g:p:r:w:a:GO:I:F:T", long_options, &option_index)) != -1){
        switch (x) {
        case 'h':
            help();
            exit(0);
            break;
        case 'c':
			globals._continually = 1;
            break;
        case 'd':
            rtl._device = atoi(optarg);
            break;
        case 'f':
            validateFreqOption(optarg);
            break;
        case 'g':
            rtl._gain = atoi(optarg);
            break;
        case 'p':
            rtl._ppm = atoi(optarg);
            break;
        case 'r':
            rtl._rate = atoi(optarg);
            break;
        case 'w':
            fftw._window_type = validateFFTWindowingOption(optarg);
            break;
        case 'a':
            fftw._avg = atoi(optarg);
            break;
		case 'G':
			globals._gnuplot = 1;
			globals._output_format = _TXT_FORMAT_;
			break;
		case 'O':
			globals._output_file = optarg;
			break;
		case 'I':
			globals._interval_seconds = atoi(optarg);
			break;
		case 'F':
			if(!globals._gnuplot){globals._output_format = validateFormatOption(optarg);}
			break;
		case 'T':
			globals._truncatefile = 1;
			break;
        default:
            help();
            exit(0);
        }
	}
	
	rtlsdr_init();
	VarAndBuffInit();
	bool cont = 1;
	while (cont) {
		time_t start_time = time(NULL);
		ReadAndRun();
		if (globals._interval_seconds > 0) {
			time_t elapsed_time = time(NULL) - start_time;
			int sleep_time = globals._interval_seconds - elapsed_time;
			if (sleep_time > 0) {
				sleep(sleep_time);
			} else {
				fprintf(stderr, "#Error: Execution time of ReadAndRun exceeds the specified interval.\n");
			}		
		}
		
		cont = globals._continually;
	}
	
	VarAndBuffDeInit();
	
	return 0;
}
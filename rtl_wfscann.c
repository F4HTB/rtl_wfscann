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
#include <sys/time.h>

//##################################main variables
struct globalsinfo{
	bool _continually = 0;
	bool _gnuplot = 0;
	bool _truncatefile = 0;
	int _output_format = _TXT_FORMAT_;
	char *_output_file = NULL;
	int _interval_seconds = -1;
	time_t _start_time;
	time_t _elapsed_time;
	bool _test = 0;
	struct timeval _test_start_time, _test_end_time;
}globals;

//##################################RTL variables
static rtlsdr_dev_t *rtl_dev;
struct rtlinfo{
	int _rate = 2048000;
	int _device = 0; 
	int _freq_start = 144000000; 
	int _freq_stop = 146000000; 
	int _freq_start_with_overlap = 144000000; 
	int _freq_stop_real = 146048000; 
	int _freq_stop_adjusted_with_overlap = 146048000;
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
	int _FFT_overlap = 0;
	int _SPS_numSlices = 0;
	int _SPS_samplesPerSlice = 0;
	fftw_complex **_all_fftin;
	fftw_complex *_fftout;
	// double **_pwr_fftw_out;
	double *_pwr_fftw_out;
	double *_dB_fftw_out;
	int _FFT_remove_binperslice = 0;
	int _FFT_final_binperslice = 0;
	int _FFT_final_binTotal  = 0;
	int _FFT_final_binTotal_Without_excess  = 0;
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
		fprintf(stderr, "\n#No supported devices found.\n");
		exit(1);
	}
	if(globals._test)
	{
		printf("#Starting rtl_map ~");
		printf("\n#Found %d device(s):", device_count);
		for(int i = 0; i < device_count; i++){
			printf(" %d: %s ", i, rtlsdr_get_device_name(i));
		}
		printf("\n");
	}
	int dev_open = rtlsdr_open(&rtl_dev, rtl._device);
	if (dev_open < 0) {
		fprintf(stderr, "\n#Failed to open RTL-SDR device #%d\n", rtl._device);
		exit(1);
	}else{
		if(globals._test){printf("#Using device: %d\n", dev_open);}
	}
	if(!rtl._gain){
		rtlsdr_set_tuner_gain_mode(rtl_dev, rtl._gain);
		if(globals._test){printf("#Gain mode set to auto.\n");}
	}else{
		rtlsdr_set_tuner_gain_mode(rtl_dev, 1);
		int gain_count = rtlsdr_get_tuner_gains(rtl_dev, NULL);
		if(globals._test){printf("#Supported gain values (%d): ", gain_count);}
		int gains[gain_count], supported_gains = rtlsdr_get_tuner_gains(rtl_dev, gains), target_gain=0;
		for (int i = 0; i < supported_gains; i++){
			if (gains[i] < rtl._gain*10){target_gain = gains[i];}
			if(globals._test){printf(" %.2f ", gains[i]/10.0);}
		}
		if(globals._test){
			printf("\n");
			printf("#Gain set to %.1f\n", target_gain / 10.0);
		}
		rtlsdr_set_tuner_gain(rtl_dev, target_gain);
	}
	rtlsdr_set_freq_correction(rtl_dev, rtl._ppm);
	rtlsdr_set_center_freq(rtl_dev, rtl._freq_start);
	rtlsdr_set_sample_rate(rtl_dev, rtl._rate);
	if(globals._test){
		printf("#Frequency start set to %d Hz.\n", rtl._freq_start);
		printf("#Sampling at %d sps\n", rtl._rate);
	}
	int r = rtlsdr_reset_buffer(rtl_dev);
	if (r < 0){
		fprintf(stderr, "\n#Failed to reset buffers.\n");
		return 1;
	}
	return 0;
}

//##################################function for images
typedef struct tagBITMAPHEADER
{
    uint16_t   bfType; // 2  /* Magic identifier */
    uint32_t   bfSize; // 4  /* File size in bytes */
    uint16_t   bfReserved1; // 2
    uint16_t   bfReserved2; // 2
    uint32_t   bfOffBits; // 4 /* Offset to image data, bytes */ 
    uint32_t    biSize; // 4 /* Header size in bytes */
    int32_t     biWidth; // 4 /* Width of image */
    int32_t     biHeight; // 4 /* Height of image */
    uint16_t    biPlanes; // 2 /* Number of colour planes */
    uint16_t    biBitCount; // 2 /* Bits per pixel */
    uint32_t    biCompress; // 4 /* Compression type */
    uint32_t    biSizeImage; // 4 /* Image size in bytes */
    int32_t     biXPelsPerMeter; // 4
    int32_t     biYPelsPerMeter; // 4 /* Pixels per meter */
    uint32_t    biClrUsed; // 4 /* Number of colours */ 
    uint32_t    biClrImportant; // 4 /* Important colours */ 
} __attribute__((packed)) BITMAPHEADER;

struct ImageData{
	int32_t _height = 1;
    int32_t _width;
    uint8_t* _data; // Tableau 1D de pixels (4 canaux par pixel)
	int _max_value_db = 0;
	int _min_value_db = -110;
	double _scaleFttToPng; //255.0 / (max_value_db - min_value_db)
} Image;

void Init_image_data(int32_t width)
{
	Image._height = 1;
	Image._width = width;
	Image._data = (uint8_t*)malloc(Image._width * 3 * sizeof(uint8_t));
	Image._scaleFttToPng = 255.0 / (Image._max_value_db - (Image._min_value_db));
}

void convert_db_data_to_image_png(FILE *fp,double* fftw_data)
{
	BITMAPHEADER bmp;
	int size = 0;
	
    for (int32_t i = 0; i < Image._width; ++i) {
        int index = (int)(((fftw_data[i] - Image._min_value_db) * Image._scaleFttToPng));
        index = (index < 0) ? 0 : ((index > 255) ? 255 : index);
        int32_t base_index = i * 3;
        Image._data[base_index] = websdrWaterfall[index][2];
        Image._data[base_index + 1] = websdrWaterfall[index][1];
        Image._data[base_index + 2] = websdrWaterfall[index][0];
    }
	
	size_t bytesRead = fread(&bmp, 1, sizeof(BITMAPHEADER), fp);
	
	if (bytesRead && bmp.bfType == 0x4D42) {
		bmp.biHeight = bmp.biHeight-1;
		size = bmp.biWidth * bmp.biHeight * 3;
        bmp.bfSize = size + sizeof(BITMAPHEADER);
		bmp.bfOffBits = sizeof(BITMAPHEADER);
		bmp.biSizeImage = size;
		
	}
	else{
		size = Image._width * Image._height * 3;
		bmp.bfType = 0x4D42; 
		bmp.bfSize= size + sizeof(BITMAPHEADER); 
		bmp.bfReserved1 = bmp.bfReserved2 = 0;
		bmp.bfOffBits = sizeof(BITMAPHEADER);
		bmp.biSize = 40;
		bmp.biWidth = Image._width;
		bmp.biHeight = -Image._height;
		bmp.biPlanes = 1;
		bmp.biBitCount = 24; 
		bmp.biCompress = 0;
		bmp.biSizeImage = size;
		bmp.biXPelsPerMeter = 0;
		bmp.biYPelsPerMeter = 0;
		bmp.biClrUsed = 0 ;
		bmp.biClrImportant = 0;
	}
	
	fseek(fp, 0, SEEK_SET);
    fwrite(&bmp, 1, sizeof(BITMAPHEADER), fp);
	fseek(fp, 0, SEEK_END);
    fwrite(Image._data, 1, bmp.biWidth * 3, fp);
}

//##################################function for options
void help() {
    printf("Usage: rtl_wfscann [OPTIONS]\n");
    printf("Options:\n");
    printf("  -h, --help               Display this help message\n");
    printf("  -c, --continuous         Run continuously\n");
    printf("  -d, --device DEVICE_ID   Set RTL device ID\n");
    printf("  -f, --freq FREQ_RANGE    Set frequency range option start:stop:fft bin width or start:stop:fft bin width:overlap in percent of sample rate (e.g., 125000000:155000000:100 or 143000000:146000000:1000:10)\n");
    printf("  -g, --gain GAIN          Set RTL gain, max 50\n");
    printf("  -p, --ppm PPM            Set RTL PPM\n");
    printf("  -r, --rate RATE          Set RTL rate\n");
    printf("  -w, --window WINDOW      Set FFT window\n");
    printf("  -a, --avg AVG            Set FFT averaging\n");
    printf("  -G, --gnuplot            Include Gnuplot option and force FORMAT to TXT (e.g., plot results)\n");
    printf("  -O, --output FILE        Set output file. With -T creat multiple file one per pass with timestamp in name file. It creat .info file with informations\n");
    printf("  -I, --interval SECONDS   Set measurement interval in seconds 0 or -1 is disable as default\n");
    printf("  -F, --format FORMAT      Set output format (TXT, IMG or IMG:db_max:db_min, RAW) default is TXT, IMG:db_max:db_min is for scale on pixel color default is 0:-110\n");
    printf("  -T, --truncatefile       Set truncated output filename like timestamp_filename.txt\n");    
	printf("\n\nExample:\n");
	
	printf("One pass of scan 143Mhz to 147Mhz 2khz bandwhith resolution with rtl key gain to 10dB:\n");
    printf("./rtl_wfscann -f 143000000:147000000:2000 -g 10\n");
	
	printf("Continuous scan 143Mhz to 147Mhz 2khz bandwhith resolution, rtl key gain to 10dB, Ouput \"Freq dBm\" format to output.txt with 5s interval:\n");
    printf("./rtl_wfscann -c -f 143000000:147000000:2000 -g 10 -F TXT -o output.txt -i 5 \n");
	
	printf("Same with creat multiple truncated file and add timestamp at file name like 1705681356_output.txt:\n");
	printf("./rtl_wfscann -c -f 143000000:147000000:2000 -g 10 -F TXT -o output.txt -i 5 -T\n");
	
	printf("One pass of scan from 143Mhz to 146Mhz, 1khz bin with, 10%% of sample rate oversample, rtl gain to 5dB\naverage mesure to 2 pass, output format TXT with gnuplot command to refresh:\n");
	printf("./rtl_wfscann -f 143000000:146000000:1000:10 -g 5 -a 2 -F TXT -G | gnuplot -p\n");
	
	printf("Continuous scan with  rtl gain to 10db, average mesure to 4 pass, bmp image output and output file to /tmp/output.bmp:\n");
	printf("./rtl_wfscann -c -f 143000000:146000000:1000:10 -g 10 -a 4 -F IMG -o /tmp/output.bmp\n");
	
	printf("Same with bmp image output and set fixed scale of colors:\n");
	printf("./rtl_wfscann -c -f 143000000:146000000:1000:10 -g 10 -a 4 -F IMG:0:-120 -o /tmp/output.bmp\n");

	printf("Same with multiple bmp file:\n");
	printf("./rtl_wfscann -c -f 143000000:146000000:1000:10 -g 10 -a 4 -F IMG:0:-120 -o /tmp/output.bmp -T\n");
	
	printf("\nNote:\n");
	printf("TXT is for view freq v value, IMG is more viewable, RAW is bether for save values.\n");
	printf("For IMG, no auto scale function because loss of performance.\n");
}

void validateFreqOption(char *freqRange) {
    char *token = strtok(freqRange, ":");

    if (token != NULL) {
        rtl._freq_start = atoi(token);
    } else {
        fprintf(stderr, "\n#Error: Missing value for --freq\n");
        exit(1);
    }

    token = strtok(NULL, ":");

    if (token != NULL) {
        rtl._freq_stop = atoi(token);
    } else {
        fprintf(stderr, "\n#Error: Missing value for --freq\n");
        exit(1);
    }

    token = strtok(NULL, ":");

    if (token != NULL) {
        fftw._bin_width = atoi(token);
    } else {
        fprintf(stderr, "\n#Error: Missing value for --freq\n");
        exit(1);
    }
	
    token = strtok(NULL, ":");

    if (token != NULL) {
        fftw._FFT_overlap = atoi(token);
    } 

	rtl._freq_start_with_overlap = rtl._freq_start - (rtl._rate*(fftw._FFT_overlap/2))/100; 

    if (rtl._freq_start_with_overlap < 0 || rtl._freq_stop < 0 || fftw._bin_width <= 0 || rtl._freq_start >= rtl._freq_stop || fftw._FFT_overlap < 0 || fftw._FFT_overlap > 25) {
        fprintf(stderr, "\n#Error: Values for --freq are invalid\n");
        exit(1);
    }

	int adjustedBandwidth = 0;
	while((rtl._freq_start+adjustedBandwidth*(rtl._rate-(rtl._rate*fftw._FFT_overlap)/100))<rtl._freq_stop)adjustedBandwidth++;
	rtl._freq_stop_adjusted_with_overlap = rtl._freq_start+adjustedBandwidth*(rtl._rate-(rtl._rate*fftw._FFT_overlap)/100);
	
	rtl._freq_stop_real = rtl._freq_stop_adjusted_with_overlap + (rtl._rate*(fftw._FFT_overlap/2))/100; 
}

void validateFormatOption(const char *format) {
    if (strcmp(format, "TXT") == 0) {
        globals._output_format = _TXT_FORMAT_;
    } else if (strcmp(format, "IMG") == 0) {
        globals._output_format = _IMG_FORMAT_;
    } else if (strncmp(format, "IMG:", 4) == 0) {
		sscanf(format, "IMG:%d:%d", &Image._max_value_db, &Image._min_value_db);
        globals._output_format = _IMG_FORMAT_;
    } else if (strcmp(format, "RAW") == 0) {
        globals._output_format = _RAW_FORMAT_;
    } else {
        fprintf(stderr, "\n#Unsupported output format: %s (TXT|IMG|RAW)\n", format);
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
        fprintf(stderr, "\n#Unsupported FFT windowing option: %s (rectangle|hanning|hamming|blackman|triangulaire|bartlett)\n", window);
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
        fprintf(stderr, "\n#Memory allocation error for FFT window.\n");
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
            fprintf(stderr, "\n#Unsupported window type.\n");
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

void ReadAndRun(){

	int n_read=0;
	for (int j = 0; j < fftw._SPS_numSlices; j=j+fftw._avg) {
		rtlsdr_set_center_freq(rtl_dev, rtl._currentFreq[j]);
		rtlsdr_reset_buffer(rtl_dev);
		for (int i = 0; i < fftw._avg; i++) {
			if (rtlsdr_read_sync(rtl_dev, rtl._samplesBuffer[j+i], fftw._SPS_samplesPerSlice, &n_read) < 0) {
				fprintf(stderr, "\n#Memory allocation error for sample buffer.\n");		
			}
			if(n_read!=fftw._SPS_samplesPerSlice){
				fprintf(stderr, "\n#Error reading samples %d/%d.\n",n_read,fftw._SPS_samplesPerSlice);		
			}
		}
	}

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
	
	memset(fftw._pwr_fftw_out, 0, fftw._FFT_final_binTotal * sizeof(double));
		
    for (int i = 0; i < numThreads; i++) {
        pthread_join(threads[i], NULL);
    }

	for (int i = 0; i < fftw._SPS_numSlices; i++) {
		fftw_execute_dft(fftPlan, fftw._all_fftin[i], fftw._fftout);
		int index = i / fftw._avg;
		int k = 0;
		for (int j = fftw._FFT_remove_binperslice/2; j < fftw._FFT_samplesPerSlice-fftw._FFT_remove_binperslice/2; j++) {
			if(j<(fftw._FFT_samplesPerSlice/2)) {k=(fftw._FFT_samplesPerSlice/2)-j;}
			else{k=(fftw._FFT_samplesPerSlice)-(j-(fftw._FFT_samplesPerSlice/2));}
			fftw._pwr_fftw_out[index*fftw._FFT_final_binperslice + (j-fftw._FFT_remove_binperslice/2)] += fftw._fftout[k][_Q_] * fftw._fftout[k][_Q_] + fftw._fftout[k][_I_] * fftw._fftout[k][_I_];
		}
	}
	
	for (int i = 0; i < fftw._FFT_final_binTotal_Without_excess; i++) {
		fftw._dB_fftw_out[i] = round(((10*log10(fftw._pwr_fftw_out[i]/fftw._avg))+ db_adc_rtl_factor)* 100) / 100;
	}
	
	if(globals._output_file) {
		char* filename = NULL;
		FILE *fichier = NULL;
		if(globals._truncatefile){
			filename = generateFileNameWithTimestamp(globals._output_file);
		}else{
			filename = globals._output_file;
		}
		switch (globals._output_format) {
			case _TXT_FORMAT_:
				fichier = fopen(filename, "w");
				for (int i = 0; i < fftw._FFT_final_binTotal_Without_excess; i++) {
					fprintf(fichier, "%d %f\n", rtl._freq_start + (i * fftw._bin_width), fftw._dB_fftw_out[i]);
				}
				fclose(fichier);
				break;
			case _IMG_FORMAT_:
				fichier = fopen(filename, "r+b");
				if (!fichier) {fichier = fopen(filename, "w+b");}
				convert_db_data_to_image_png(fichier,fftw._dB_fftw_out);
				fclose(fichier);
				break;
			case _RAW_FORMAT_:
				fichier = fopen(filename, "wb");
				fwrite(fftw._dB_fftw_out, sizeof(double), fftw._FFT_final_binTotal_Without_excess, fichier);
				fclose(fichier);
				break;
			default:
				break;
			}
	}
	else{
		switch (globals._output_format) {
			case _TXT_FORMAT_:
				if(globals._gnuplot){printf("e\nplot \"-\"\n");}
				for (int i = 0; i < fftw._FFT_final_binTotal_Without_excess; i++) {
					printf("%d %f\n", rtl._freq_start + (i * fftw._bin_width), fftw._dB_fftw_out[i]); 
				}
				break;
			case _RAW_FORMAT_:
				fwrite(fftw._dB_fftw_out, sizeof(double), fftw._FFT_final_binTotal_Without_excess, stdout);
				break;
			default:
				break;
			}
	}
}

void VarAndBuffInit(){
	
	fftPlan = fftw_plan_dft_1d((rtl._rate / fftw._bin_width), NULL, NULL, FFTW_BACKWARD, FFTW_ESTIMATE);

	fftw._FFT_numSlices = (rtl._freq_stop_adjusted_with_overlap - rtl._freq_start) / (rtl._rate-(rtl._rate*fftw._FFT_overlap)/100);
	fftw._FFT_samplesPerSlice = (rtl._rate / fftw._bin_width); // /2 because fftw_complex is I/Q
	fftw._SPS_numSlices = fftw._FFT_numSlices * fftw._avg;
	fftw._SPS_samplesPerSlice = fftw._FFT_samplesPerSlice * 2;// *2 beacause 1 sample = 2 byte I/Q
	fftw._FFT_remove_binperslice = (fftw._FFT_samplesPerSlice*fftw._FFT_overlap)/100;
	fftw._FFT_final_binperslice = fftw._FFT_samplesPerSlice - fftw._FFT_remove_binperslice;
	fftw._FFT_final_binTotal = fftw._FFT_numSlices*fftw._FFT_final_binperslice;
	fftw._FFT_final_binTotal_Without_excess = (fftw._FFT_final_binTotal - (rtl._freq_stop_adjusted_with_overlap - rtl._freq_stop)/ fftw._bin_width);

	if(globals._test)
	{
		printf("\n");
		printf("##################################################################################\n");
		printf("#\n");
		printf("#1 Frequency start wanted : %dhz\n", rtl._freq_start);
		printf("#2 Frequency start adjusted with overlap = real = \"Frequency start wanted\" - ((\"sample rate\"*(\"fft overlappercent\"/2))/100)\n");
		printf("#2 %d - ((%d*(%d/2))/100) = %dhz\n",rtl._freq_start, rtl._rate, fftw._FFT_overlap, rtl._freq_start_with_overlap);
		printf("#3 Frequency stop wanted : %dhz\n", rtl._freq_stop);
		
		printf("#4 Frequency stop adjusted with overlap = \"Frequency start wanted\" + \"number of slice\" * \"Frequency bandwidth per slice\"\n");	
		printf("#4 Frequency stop adjusted with overlap = %d + %d * %d = %dhz\n",rtl._freq_start,fftw._FFT_numSlices,(rtl._rate-(rtl._rate*fftw._FFT_overlap)/100),rtl._freq_stop_adjusted_with_overlap);	
		
		printf("#5 Frequency stop real = \"Frequency start adjusted with overlap\" + \"number of slice\" * \"sample rate\" - (\"number of slice\" - 1) * \"frequency overlaped\" \n");
		printf("#5 Frequency stop real = %d + %d * %d - %d * %d = %dhz\n",rtl._freq_start_with_overlap,fftw._FFT_numSlices,rtl._rate,fftw._FFT_numSlices-1,((rtl._rate*fftw._FFT_overlap)/100),rtl._freq_stop_real);
		
		printf("#5 Frequency overlap per scan bandwith (rlt sample rate): %d%% %dhz\n",fftw._FFT_overlap,(rtl._rate*(fftw._FFT_overlap))/100);
		printf("#\n");
		printf("#2  1                                  3    4   \n");
		printf("#|__|---------|____|---------|____|----|----|__|\n");
		printf("# 5/           5/\\5           5/\\5           \\5 \n");
		printf("#\n");
		printf("#RTL number of scan:%d with avg: %d\n",fftw._FFT_numSlices,fftw._SPS_numSlices);
		printf("#RTL number of sample per scan:%d with avg: %d\n",fftw._SPS_samplesPerSlice/2, (fftw._SPS_samplesPerSlice/2) * fftw._avg);
		printf("#RTL total number of sample %d\n",fftw._SPS_numSlices * fftw._avg * fftw._SPS_samplesPerSlice/2);
		printf("#\n");
		printf("#FFT number of fftw to compute %d with avg %d\n",fftw._FFT_numSlices, fftw._FFT_numSlices * fftw._avg);
		printf("#FFT number of bin per fftw %d with avg %d\n",fftw._FFT_samplesPerSlice, fftw._FFT_samplesPerSlice * fftw._avg);	
		printf("#FFT total number of bin with avg %d\n", fftw._FFT_samplesPerSlice * fftw._avg * fftw._FFT_numSlices);	
		printf("#FFT number of bin to remove per fft for overlap %d\n",fftw._FFT_remove_binperslice);
		printf("#FFT number of fft bin after overlap %d\n",fftw._FFT_final_binperslice);	
		printf("#FFT total number of fft bin after overlap %d\n",fftw._FFT_final_binTotal);	
		printf("#FFT total number of fft bin to stop frequency %d\n",fftw._FFT_final_binTotal_Without_excess);
		printf("#FFT total number of fft bin that not finaly computed %d\n",fftw._FFT_final_binTotal_Without_excess-fftw._FFT_final_binTotal);
		printf("#\n");		
		printf("##################################################################################\n\n");
	}
	
	rtl._currentFreq = (int *)malloc(fftw._SPS_numSlices * sizeof(int));
	if (rtl._currentFreq == NULL) {
		fprintf(stderr, "\n#Memory allocation error for frequencies index.\n");
		exit(1);
	}
	for (int j = 0; j < fftw._SPS_numSlices; j=j+fftw._avg) {
		rtl._currentFreq[j] = rtl._freq_start + (j/fftw._avg * (rtl._rate - (rtl._rate*fftw._FFT_overlap)/100)) + ((rtl._rate - (rtl._rate*fftw._FFT_overlap)/100) / 2);
	}
	
	fftw._window = generateFFTWindow(fftw._FFT_samplesPerSlice,fftw._window_type);	

	rtl._samplesBuffer = (uint8_t **)malloc(fftw._SPS_numSlices * sizeof(uint8_t *));
    if (rtl._samplesBuffer == NULL) {
        fprintf(stderr, "\n#Memory allocation error for sample buffer.\n");
        exit(1);
    }
	
	for (int i = 0; i < fftw._SPS_numSlices; i++) {
		rtl._samplesBuffer[i] = (uint8_t *)malloc(fftw._SPS_samplesPerSlice * sizeof(uint8_t));
        if (rtl._samplesBuffer[i] == NULL) {
            fprintf(stderr, "\n#Memory allocation error for sample buffer.\n");
            exit(1);
        }
	}
	
	fftw._all_fftin = (fftw_complex **)fftw_malloc(fftw._SPS_numSlices * sizeof(fftw_complex *));
	if (fftw._all_fftin == NULL) {
        fprintf(stderr, "\n#Memory allocation error for fft input buffer.\n");
        exit(1);
    }
	for (int i = 0; i < fftw._SPS_numSlices; i++) {
		fftw._all_fftin[i]= (fftw_complex *)fftw_malloc(fftw._FFT_samplesPerSlice * sizeof(fftw_complex));
        if (fftw._all_fftin[i] == NULL) {
            fprintf(stderr, "\n#Memory allocation error for fft input buffer.\n");
            exit(1);
        }
	}
	
	fftw._fftout = (fftw_complex *)fftw_malloc(fftw._FFT_samplesPerSlice * sizeof(fftw_complex));
    if (fftw._fftout == NULL) {
        fprintf(stderr, "\n#Memory allocation error for fft output buffer\n");
        exit(1);
    }
	
	fftw._pwr_fftw_out = (double *)calloc(fftw._FFT_final_binTotal,sizeof(double));
	if (fftw._pwr_fftw_out == NULL) {
		fprintf(stderr, "\n#Memory allocation error for fft power values.\n");
		exit(1);
	}
	
	numThreads = sysconf(_SC_NPROCESSORS_ONLN);
	threads = (pthread_t*)malloc(numThreads * sizeof(pthread_t));
	threadArgsArray = (struct ThreadArgs*)malloc(numThreads * sizeof(struct ThreadArgs));
	if (threadArgsArray == NULL) {
		fprintf(stderr, "\n#Thread allocation error for prep threads.\n");
		exit(1);
	}
	
	fftw._dB_fftw_out = (double *)calloc(fftw._FFT_final_binTotal_Without_excess,sizeof(double));
	if (fftw._dB_fftw_out == NULL) {
		fprintf(stderr, "\n#Memory allocation error for fft power values in dB.\n");
		exit(1);
	}
	
	if(globals._output_file) {
		
		char* new_filename = (char*)malloc(strlen(globals._output_file) + strlen(".info") + 1);
		strcpy(new_filename, globals._output_file);strcat(new_filename, ".info");
		
		FILE *fichier = fopen(new_filename, "w");
		printf("#File name: %s\n", globals._output_file);
		if (fichier == NULL) {
			fprintf(stderr, "\n#Impossible to open the file.\n");
			exit(1);
		}
		fprintf(fichier, "freq_start:%d\nfreq_stop:%d\nbin_bandwith:%d\ntotal_samples:%d\n", rtl._freq_start, rtl._freq_stop,fftw._bin_width,fftw._FFT_final_binTotal_Without_excess);
		fclose(fichier);
		
		if(globals._output_format == _IMG_FORMAT_){
				Init_image_data(fftw._FFT_final_binTotal_Without_excess);
		}
	}
	

}

void VarAndBuffDeInit(){
	
	free(threadArgsArray);
	free(threads);
	
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
			{"test", no_argument, 0, 't'},
			{0, 0, 0, 0}
        };
	int option_index = 0;

    while ((x = getopt_long(argc, argv, "hcd:f:g:p:r:w:a:GO:I:F:Tt", long_options, &option_index)) != -1){
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
			if(!globals._gnuplot){validateFormatOption(optarg);}
			break;
		case 'T':
			globals._truncatefile = 1;
			break;
		case 't':
			globals._test = 1;
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
		
		if(globals._test){
			gettimeofday(&globals._test_start_time, NULL);
		}
		
		globals._start_time = time(NULL);
		ReadAndRun();
		if (globals._interval_seconds > 0) {
			globals._elapsed_time = time(NULL) - globals._start_time;
			int sleep_time = globals._interval_seconds - globals._elapsed_time;
			if (sleep_time >= 0) {
				sleep(sleep_time);
			} else {
				fprintf(stderr, "\n#Error: Execution time of pass exceeds the specified %ds interval: %ds .\n",globals._interval_seconds, globals._interval_seconds-sleep_time);
			}		
		}
		
		cont = globals._continually;
		
		if(globals._test){
			gettimeofday(&globals._test_end_time, NULL);
			long elapsed_time = (globals._test_end_time.tv_sec - globals._test_start_time.tv_sec) * 1000000L + (globals._test_end_time.tv_usec - globals._test_start_time.tv_usec);
			double elapsed_time_seconds = (double)elapsed_time / 1000000.0;
			printf("#Temps d'exécution : %.6f microsecondes\n", elapsed_time_seconds);
		}

	}
	
	VarAndBuffDeInit();
	
	return 0;
}
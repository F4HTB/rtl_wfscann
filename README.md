# rtl_wfscann
SDR_POWER like suite of tools command line and with web interface.
```
Usage: rtl_wfscann [OPTIONS]
Options:
  -h, --help               Display this help message
  -c, --continuous         Run continuously
  -d, --device DEVICE_ID   Set RTL device ID
  -f, --freq FREQ_RANGE    Set frequency range option start:stop:fft bin width or start:stop:fft bin width:overlap in percent of sample rate (e.g., 125000000:155000000:100 or 143000000:146000000:1000:10)
  -g, --gain GAIN          Set RTL gain, max 50
  -p, --ppm PPM            Set RTL PPM
  -r, --rate RATE          Set RTL rate
  -w, --window WINDOW      Set FFT window
  -a, --avg AVG            Set FFT averaging
  -G, --gnuplot            Include Gnuplot option and force FORMAT to TXT (e.g., plot results)
  -O, --output FILE        Set output file. With -T creat multiple file one per pass with timestamp in name file. It creat .info file with informations
  -I, --interval SECONDS   Set measurement interval in seconds 0 or -1 is disable as default
  -F, --format FORMAT      Set output format (TXT, IMG or IMG:db_max:db_min, RAW) default is TXT, IMG:db_max:db_min is for scale on pixel color default is 0:-110
  -T, --truncatefile       Set truncated output filename like timestamp_filename.txt


Example:
One pass of scan 143Mhz to 147Mhz 2khz bandwhith resolution with rtl key gain to 10dB:
./rtl_wfscann -f 143000000:147000000:2000 -g 10
Continuous scan 143Mhz to 147Mhz 2khz bandwhith resolution, rtl key gain to 10dB, Ouput "Freq dBm" format to output.txt with 5s interval:
./rtl_wfscann -c -f 143000000:147000000:2000 -g 10 -F TXT -o output.txt -i 5
Same with creat multiple truncated file and add timestamp at file name like 1705681356_output.txt:
./rtl_wfscann -c -f 143000000:147000000:2000 -g 10 -F TXT -o output.txt -i 5 -T
One pass of scan from 143Mhz to 146Mhz, 1khz bin with, 10% of sample rate oversample, rtl gain to 5dB
average mesure to 2 pass, output format TXT with gnuplot command to refresh:
./rtl_wfscann -f 143000000:146000000:1000:10 -g 5 -a 2 -F TXT -G | gnuplot -p
Continuous scan with  rtl gain to 10db, average mesure to 4 pass, bmp image output and output file to /tmp/output.bmp:
./rtl_wfscann -c -f 143000000:146000000:1000:10 -g 10 -a 4 -F IMG -o /tmp/output.bmp
Same with bmp image output and set fixed scale of colors:
./rtl_wfscann -c -f 143000000:146000000:1000:10 -g 10 -a 4 -F IMG:0:-120 -o /tmp/output.bmp
Same with multiple bmp file:
./rtl_wfscann -c -f 143000000:146000000:1000:10 -g 10 -a 4 -F IMG:0:-120 -o /tmp/output.bmp -T

Note:
TXT is for view freq v value, IMG is more viewable, RAW is bether for save values.
For IMG, no auto scale function because loss of performance.

```

In progress.
The waterfall generation tool is working at this time.

The python part for the web server is in progress.
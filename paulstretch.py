#
# Paul's Extreme Sound Stretch (Paulstretch) - Python version
#
# by Nasca Octavian PAUL, Targu Mures, Romania
# http://www.paulnasca.com/
#
# http://hypermammut.sourceforge.net/paulstretch/
#
#
# Modified by Evan Viera
# this file is released under Public Domain
#

import sys
import scipy.io.wavfile
import wave
from numpy import *
from optparse import OptionParser
import os

# Function to load WAV file
def load_wav(filename):
    try:
        wavedata = scipy.io.wavfile.read(filename)
        samplerate = int(wavedata[0])
        smp = wavedata[1] * (1.0/32768.0)
        smp = smp.transpose()
        if len(smp.shape) == 1:  # Convert to stereo if it's not
            smp = tile(smp, (2, 1))
        return (samplerate, smp)
    except:
        print ("Error loading wav: "+filename)
        return None

# Function to optimize window size
def optimize_windowsize(n):
    orig_n = n
    while True:
        n = orig_n
        # Reduce n by its prime factors 2, 3, and 5
        while (n%2) == 0:
            n /= 2
        while (n%3) == 0:
            n /= 3
        while (n%5) == 0:
            n /= 5
        if n < 2:
            break
        orig_n += 1
    return orig_n

# Function to perform Paulstretch
def paulstretch(samplerate, smp, stretch, windowsize_seconds, outfilename):
    nchannels = smp.shape[0]
    outfile = wave.open(outfilename, "wb")
    outfile.setsampwidth(2)
    outfile.setframerate(samplerate)
    outfile.setnchannels(nchannels)

    # Initialize window size
    windowsize = int(windowsize_seconds * samplerate)
    if windowsize < 16:
        windowsize = 16
    windowsize = optimize_windowsize(windowsize)
    windowsize = int(windowsize / 2) * 2
    half_windowsize = int(windowsize / 2)

    # Prepare the sample
    nsamples = smp.shape[1]
    end_size = int(samplerate * 0.05)
    if end_size < 16:
        end_size = 16
    smp[:, nsamples-end_size:nsamples] *= linspace(1, 0, end_size)

    # Initialize position variables
    start_pos = 0.0
    displace_pos = (windowsize * 0.5) / stretch

    # Create window
    window = pow(1.0 - pow(linspace(-1.0, 1.0, windowsize), 2.0), 1.25)
    old_windowed_buf = zeros((2, windowsize))

    while True:
        istart_pos = int(floor(start_pos))
        buf = smp[:, istart_pos:istart_pos + windowsize]
        if buf.shape[1] < windowsize:
            buf = append(buf, zeros((2, windowsize-buf.shape[1])), 1)
        buf = buf * window

        # Get the amplitudes of the frequency components and discard the phases
        freqs = abs(fft.rfft(buf))

        # Randomize phases
        ph = random.uniform(0, 2*pi, (nchannels, freqs.shape[1])) * 1j
        freqs = freqs * exp(ph)

        buf = fft.irfft(freqs)
        buf *= window

        # Output audio
        output = buf[:, 0:half_windowsize] + old_windowed_buf[:, half_windowsize:windowsize]
        old_windowed_buf = buf
        output[output > 1.0] = 1.0
        output[output < -1.0] = -1.0
        outfile.writeframes(int16(output.ravel('F') * 32767.0).tobytes())

        start_pos += displace_pos
        if start_pos >= nsamples:
            print ("100 %")
            break
        print ("%d %% \r" % int(100.0 * start_pos / nsamples), end='', flush=True)

    outfile.close()

# Start of the main program
print ("Paul's Extreme Sound Stretch (Paulstretch) - Python version 20141220")
print ("by Nasca Octavian PAUL, Targu Mures, Romania\n")

# Create command line option parser
parser = OptionParser(usage="usage: %prog [options] input_wav")
parser.add_option("-s", "--stretch", dest="stretch",help="stretch amount (1.0 = no stretch)",type="float",default=8.0)
parser.add_option("-w", "--window_size", dest="window_size",help="window size (seconds)",type="float",default=0.25)
(options, args) = parser.parse_args()

# Check for valid input arguments
if (len(args) < 1) or (options.stretch <= 0.0) or (options.window_size <= 0.001):
    print ("Error in command line parameters. Run this program with --help for help.")
    sys.exit(1)

print ("stretch amount = %g" % options.stretch)
print ("window size = %g seconds" % options.window_size)

# Load WAV file
(samplerate, smp) = load_wav(args[0])

# Generate output file name based on input name and parameters
input_name, input_ext = os.path.splitext(args[0])
stretch_str = str(options.stretch).replace('.', 'p')
window_size_str = str(options.window_size).replace('.', 'p')
output_name = f"{input_name}_stretched_s{stretch_str}_w{window_size_str}{input_ext}"

# Run Paulstretch
paulstretch(samplerate, smp, options.stretch, options.window_size, output_name)

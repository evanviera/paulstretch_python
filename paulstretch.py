#
# Paul's Extreme Sound Stretch (Paulstretch) - Python version
#
# by Nasca Octavian PAUL, Targu Mures, Romania
# http://www.paulnasca.com/
#
# http://hypermammut.sourceforge.net/paulstretch/
#
#
# modified by Evan Viera 
#
# this file is released under Public Domain
#

import sys
import scipy.io.wavfile
import wave
from numpy import *
from optparse import OptionParser

# Function to load a wav file
def load_wav(filename):
    try:
        # Read the wav file
        wavedata=scipy.io.wavfile.read(filename)
        # Get the sample rate
        samplerate=int(wavedata[0])
        # Normalize the samples
        smp=wavedata[1]*(1.0/32768.0)
        smp=smp.transpose()
        # If mono, convert to stereo
        if len(smp.shape)==1:
            smp=tile(smp,(2,1))
        return (samplerate,smp)
    except:
        print ("Error loading wav: "+filename)
        return None


# Function to optimize the window size
def optimize_windowsize(n):
    orig_n=n
    while True:
        n=orig_n
        # Divide by 2 as long as n is even
        while (n%2)==0:
            n/=2
        # Divide by 3 as long as n is divisible by 3
        while (n%3)==0:
            n/=3
        # Divide by 5 as long as n is divisible by 5
        while (n%5)==0:
            n/=5

        if n<2:
            break
        orig_n+=1
    return orig_n


# The main function for sound stretching
def paulstretch(samplerate,smp,stretch,windowsize_seconds,outfilename):
    # Number of channels (stereo)
    nchannels=smp.shape[0]

    # Initialize the output file
    outfile=wave.open(outfilename,"wb")
    outfile.setsampwidth(2)
    outfile.setframerate(samplerate)
    outfile.setnchannels(nchannels)

    # Adjust the window size
    windowsize=int(windowsize_seconds*samplerate)
    if windowsize<16:
        windowsize=16
    windowsize=optimize_windowsize(windowsize)
    windowsize=int(windowsize/2)*2
    half_windowsize=int(windowsize/2)

    # Correct the end of the smp
    nsamples=smp.shape[1]
    end_size=int(samplerate*0.05)
    if end_size<16:
        end_size=16

    smp[:,nsamples-end_size:nsamples]*=linspace(1,0,end_size)

    # Compute the displacement inside the input file
    start_pos=0.0
    displace_pos=(windowsize*0.5)/stretch

    # Create window function
    window=pow(1.0-pow(linspace(-1.0,1.0,windowsize),2.0),1.25)

    old_windowed_buf=zeros((2,windowsize))

    while True:
        # Get the windowed buffer
        istart_pos=int(floor(start_pos))
        buf=smp[:,istart_pos:istart_pos+windowsize]
        if buf.shape[1]<windowsize:
            buf=append(buf,zeros((2,windowsize-buf.shape[1])),1)
        buf=buf*window

        # Compute frequency components
        freqs=abs(fft.rfft(buf))

        # Randomize the phases
        ph=random.uniform(0,2*pi,(nchannels,freqs.shape[1]))*1j
        freqs=freqs*exp(ph)

        # Inverse FFT
        buf=fft.irfft(freqs)

        # Window the output buffer
        buf*=window

        # Overlap-add the output
        output=buf[:,0:half_windowsize]+old_windowed_buf[:,half_windowsize:windowsize]
        old_windowed_buf=buf

        # Clamp the values to -1..1 
        output[output>1.0]=1.0
        output[output<-1.0]=-1.0

        # Write the output to wav file
        outfile.writeframes(int16(output.ravel('F')*32767.0).tobytes())

        start_pos+=displace_pos
        if start_pos>=nsamples:
            print ("100 %")
            break
        print ("%d %% \r" % int(100.0*start_pos/nsamples), end='', flush=True)

    outfile.close()

print ("Paul's Extreme Sound Stretch (Paulstretch) - Python version 20141220")
print ("by Nasca Octavian PAUL, Targu Mures, Romania\n")
parser = OptionParser(usage="usage: %prog [options] input_wav output_wav")
parser.add_option("-s", "--stretch", dest="stretch",help="stretch amount (1.0 = no stretch)",type="float",default=8.0)
parser.add_option("-w", "--window_size", dest="window_size",help="window size (seconds)",type="float",default=0.25)
(options, args) = parser.parse_args()

# Check command line parameters
if (len(args)<2) or (options.stretch<=0.0) or (options.window_size<=0.001):
    print ("Error in command line parameters. Run this program with --help for help.")
    sys.exit(1)

print ("stretch amount = %g" % options.stretch)
print ("window size = %g seconds" % options.window_size)
(samplerate,smp)=load_wav(args[0])

# Run the Paulstretch algorithm
paulstretch(samplerate,smp,options.stretch,options.window_size,args[1])

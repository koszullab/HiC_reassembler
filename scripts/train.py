# Script used to train the machine learning models which are used to detect SV breakpoints.

from hiscram.detector.matrixdetector import Matrixdetector
if __name__ == "__main__":
    print("TRAIN MATRIXDETECTOR:")
    MatDetect = Matrixdetector() # Train model which detects SV breakpoints on the Hi-C matrix.
    MatDetect.train()
    MatDetect.plot()
    MatDetect.save()
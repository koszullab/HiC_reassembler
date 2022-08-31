import numpy as np
import sys

# Tool script to concatenate several datasets. For instance, if you want to establish a dataset
# based on several chromosomes, you can concatenate the datasets corresponding to each chromosome
# with this script.

imgs = []
labels = []
for i in range(1, len(sys.argv)):
    try:
        a = np.load(sys.argv[i] + "/dataset/imgs.npy")
        for img in a:
            imgs.append(img)
        b = np.load(sys.argv[i] + "/dataset/labels.npy")
    except:
        pass

imgs = np.array(imgs)
labels = np.array(labels)

print("shape of imgs = ", imgs.shape)
print("shape of labels = ", labels.shape)
np.save("compiled_set/imgs.npy", imgs)
np.save("compiled_set/labels.npy", labels)

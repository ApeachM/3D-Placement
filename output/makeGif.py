import imageio
import os

dir = 'images/'
file_type = r'bmp'
speed_spec = {'duration': 0.05}

filenames = os.listdir(dir)
filenames = sorted(filenames)
images = []

for filename in filenames:
    if filename.endswith('.{}'.format(file_type)):
        file_path = os.path.join(dir, filename)
        if "die" in file_path or "top" in file_path or "bottom" in file_path:
            continue
        images.append(imageio.imread(file_path))
        print(filename)
imageio.mimsave('{}{}.gif'.format(dir, "ani"), images, **speed_spec)

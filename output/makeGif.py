import imageio
import os

dir = 'images/'
file_type = r'bmp'
speed_spec = {'duration': 0.05}

filenames = os.listdir(dir)
filenames = sorted(filenames)
top_die_images = []
bottom_die_images = []
pseudo_die_images = []

for filename in filenames:
    if filename.endswith('.{}'.format(file_type)):
        file_path = os.path.join(dir, filename)
        if "die" in file_path or "after_hybrid_bond_generation" in file_path:
            continue
        if "top" in file_path:
            top_die_images.append(imageio.imread(file_path))
        elif "bottom" in file_path:
            bottom_die_images.append(imageio.imread(file_path))
        elif "pseudo" in file_path:
            pseudo_die_images.append(imageio.imread(file_path))
        # print(filename)
if len(top_die_images) > 1:
    imageio.mimsave('{}{}.gif'.format(dir, "top"), top_die_images, **speed_spec)
if len(bottom_die_images) > 1:
    imageio.mimsave('{}{}.gif'.format(dir, "bottom"), bottom_die_images, **speed_spec)
if len(pseudo_die_images) > 1:
    imageio.mimsave('{}{}.gif'.format(dir, "pseudo"), pseudo_die_images, **speed_spec)

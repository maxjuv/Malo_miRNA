from configuration import *
import matplotlib

def get_color(proportion , cmap = 'viridis'):

    color_map =matplotlib.cm.get_cmap(cmap)
    c = color_map(proportion)
    return c


if __name__ == '__main__':
    print(get_color(.3))

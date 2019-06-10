import os
import matplotlib.pyplot as plt

from PIL import Image


def image_to_bytes(file_name):
    cwd = os.getcwd()
    file_path = os.path.join(cwd, '..', 'files', 'images', file_name) + '.png'
    image = Image.open(file_path)

    data = list(image.getdata(band=0))  # The image is in grayscale

    byte_path = os.path.join(cwd, '..', 'files', 'bytes', file_name) + '.bt'
    try:
        byte_file = open(byte_path, 'xb')
    except FileExistsError:
        byte_file = open(byte_path, 'wb')

    for x in data:
        byte_file.write(x.to_bytes(1, 'big'))

    byte_file.close()


def bytes_to_image(file_name, size):
    cwd = os.getcwd()
    byte_path = os.path.join(cwd, '..', 'files', 'bytes', file_name) + '.bt'
    byte_file = open(byte_path, 'rb')

    image = Image.frombytes('L', size, byte_file.read(), "raw")
    file_path = os.path.join(cwd, '..', 'files', 'images', file_name) + '.png'
    try:
        file = open(file_path, 'xb')
    except FileExistsError:
        file = open(file_path, 'wb')

    image.save(file)
    file.close()
    image.close()


def convert_all(file_name, size):
    bytes_to_image(file_name + '_n', size)
    bytes_to_image(file_name + '_d', size)


def graph_ber():
    cwd = os.getcwd()
    t_path = os.path.join(cwd, '../files/graph/turbo.grp')
    l_path = os.path.join(cwd, '../files/graph/ldpc.grp')
    b_path = os.path.join(cwd, '../files/graph/base.grp')

    t_file = open(t_path, 'r')
    l_file = open(l_path, 'r')
    b_file = open(b_path, 'r')

    t_snr, t_err = [], []
    l_snr, l_err = [], []
    b_snr, b_err = [], []

    for line in t_file:
        data = line.strip('\n').split(' ')
        t_snr.append(float(data[0]))
        t_err.append(100*float(data[1]))

    for line in l_file:
        data = line.strip('\n').split(' ')
        l_snr.append(float(data[0]))
        l_err.append(100*float(data[1]))

    for line in b_file:
        data = line.strip('\n').split(' ')
        b_snr.append(float(data[0]))
        b_err.append(100*float(data[1]))

    t_file.close()
    l_file.close()
    b_file.close()

    plt.plot(t_snr, t_err, label='Turbocode')
    plt.plot(l_snr, l_err, label='Code LDPC (Hard Decoding)')
    plt.plot(b_snr, b_err, label='Sans codage')

    plt.xlabel("$\\frac{E_b}{N_0}$ (dB)")
    plt.ylabel("Taux d'erreur (%)")

    plt.legend()
    plt.show()

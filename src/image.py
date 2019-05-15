import os

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

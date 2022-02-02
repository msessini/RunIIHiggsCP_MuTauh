import os
import errno

def create_dir(file_path):
    if '/' not in file_path:
        return
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
        try:
            os.makedirs(directory)
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(directory):
                pass
            else:
                print("Bad name for directory:", file_path)
                raise


def save(canvas, name, extension='.pdf'):
    name = name.replace(' ', '').replace('&&', '')
    create_dir(name)
    canvas.SaveAs(name + extension)

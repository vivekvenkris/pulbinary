import sys
import os
import os.path
from errno import ENOENT
import subprocess
import time

def print_in_stderr(args, **kwargs):
    print(args, file=sys.stderr, **kwargs)


def convert_to_float(value):
    try:
        return float(value)
    except ValueError as e:
        print_in_stderr(e)
        return 0.0


def convert_to_int(value):
    try:
        return int(value)
    except ValueError as e:
        print_in_stderr(e)
        return 0.0


def ensure_file_exists(file_name):
    if not os.path.isfile(file_name):
        raise IOError(ENOENT, 'No file found', file_name)



def ensure_directory_exists(directory_name):
    if not os.path.isdir(directory_name):
        raise IOError(ENOENT, 'No directory found', directory_name)


def ensure_directory_is_not_empty(directory_name):
    if len(os.listdir(directory_name)) == 0:
        raise IOError(ENOENT, 'No files in directory', directory_name)

def is_number(s):
    try:
        if s is None:
            return False
        float(s)
        return True
    except ValueError:
        return False


class SubProcessRunner(object):

    def __init__(self, command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env=os.environ.copy()):
        self.stdout = stdout
        self.stderr = stderr
        self.command = command
        self.proc = None
        self.env=env

    def run(self):
        command_chunks = self.command.split()

        self.proc = subprocess.Popen(
            command_chunks, stdout=self.stdout, stderr=self.stderr, env=self.env)
        # Wait for the process to complete
        self.stdout, self.stderr = self.proc.communicate()

        if self.proc.returncode:
            print("\n RUNNING COMMAND FAILED: {}\n".format(self.command))
            print("stdout:",self.stdout)
            print("stderr:", self.stderr)






def run_process(command, env=os.environ.copy()):
    if not env:
        env = os.environ.copy()
    subprocess_runner = SubProcessRunner(command, env=env)
    subprocess_runner.run()

    return subprocess_runner.stdout.decode()

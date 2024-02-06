# Python N-body class and related files

Here's how I recommend setting up your system to use these files

* You probably have a file in your home directory called .bashrc (the dot hides the file, you might need to do "ls -a" to see it)
* Inside this file add the following lines, where "path/to/nbody" is the path from your home directory to this directory in your system (because you're using git from the command line, right?):

export PATH=/home/path/to/nbody:$PATH

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/path/to/nbody

export PYTHONPATH=$PYTHONPATH:/home/path/to/nbody

* you can then run any exectuable file you put in this directory (e.g., .py files if you make them exectuable).

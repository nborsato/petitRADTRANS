Installation
============

Prerequisites
_____________

- A working Python installation, including the numpy package.
- A Fortran compiler (e.g. “gfortran”)


Download Code
_____________

Download petitRADTRANS from `Gitlab <https://gitlab.com/mauricemolli/petitRADTRANS.git>`_, or clone it from GitLab via

.. code-block:: bash
		
   git clone git@gitlab.com:mauricemolli/petitRADTRANS.git

.. note::
   We are working on making petitRADTRANS available via pip install and anaconda, and providing as setup script. In the meantime, please carry out the (few) steps below.

Download opacity data
_____________________

.. important::
   **Dropbox temporarily bans pRT links sometimes due to bandwidth
   overuse.** If you cannot access our files, please contact us `here
   <mailto:molliere@mpia.de>`_. We are working on a better solution.

Download the `opacity and input data
<https://www.dropbox.com/s/tz6j0rowpdbs509/input_data_r.zip?dl=0>`_
(6.3 GB), unzip them, and put the "input_data" folder into the
"petitRADTRANS" folder (i.e. the same folder where the source is, if
you clone from gitlab, this should be the
petitRADTRANS folder *in* the petitRADTRANS folder). This contains the
necessary files to run petitRADTRANS, and the low resolution
(:math:`\lambda/\Delta\lambda=1000`) opacity files. The high
resolution (:math:`\lambda/\Delta\lambda=10^6`) opacity data (about
240 GB if you want to get all species) can be
accessed and downloaded `via Dropbox here`_. To
install them, create a folder called "line_by_line" in the
"input_data/opacities/lines" folder. Then put the folder of the absorber
species you downloaded in there.

.. _`via Dropbox here`: https://www.dropbox.com/sh/w7sa20v8qp19b4d/AABKF0GsjghsYLJMUJXDgrHma?dl=0
Installation
____________

- In the terminal, enter the petitRADTRANS folder containing the source
  (.py and .f90 files) (``cd petitRADTRANS``).
- Type the following in the terminal ``chmod +x make.sh``, and press Enter.
- Type the following in the terminal ``./make.sh``, and press Enter. A lot of text will appear while the Fortran subroutines are being built. If you use Anaconda, see the first installation tip below before carrying out this step.
- Type ``ls`` in the terminal and press Enter. If everything has worked you should see three files with the “.so” extension in the folder. If you are experiencing problems see the installation tips below.
- Open the “.bash_profile” or “.bashrc” file (depending on your operating system) in your home directory. Add the following as a new line to the file (you may have to use sudo for modifying this file).

.. code-block:: bash
		
   export PYTHONPATH=Path of the folder containing the petitRADTRANS source/:$PYTHONPATH

.. attention::
   Don’t forget to adapt the path in the line above :) ! If you clone petitRADTRANS from gitlab, this
   should be the path of the top-level petitRADTRANS folder
   *containing* the petitRADTRANS source folder. If you are
   uncertain what the absolute path of the folder containing the
   petitRADTRANS source folder is, then switch to that folder in the
   terminal, type “pwd”, and press Enter. You can then just copy-paste
   that path. Don’t forget to put a dash “/“
   at the end of the path.
   Close and reopen the terminal such that it will set the Python path correctly.

Installation tips
_________________

- When running “make.sh”, the shell script will call a program named “f2py”. f2py will compile the Fortran subroutines used by petitRADTRANS, and wrap them for use in Python. f2py is part of the Python numpy package. If you are running Anaconda and use different Python environments, it is important that you activate the Python environment that you want to run petitRADTRANS in before running “make.sh”. Only then will the Fortran subroutines be build with the correct f2py version (using the Python and numpy version valid within the Anaconda environment).
- If you put the petitRADTRANS folder in a protected area of your file system (like “/Applications” on Mac OS) you may have to carry out all installation steps with the “sudo” command preceding the actual terminal commands ``chmod …`` and ``./make.sh``, otherwise you may run into “Permission denied” errors. You will have to enter your password when doing this.
- Sometimes ``./make.sh`` will not work. In this case copy the lines contained within “make.sh” individually to the clipboard, paste them to the terminal, and press Enter.

Testing the installation
________________________

Open a new terminal window (this will source the ``PYTHONPATH``). Then open python and type

.. code-block:: python
		
   from petitRADTRANS import Radtrans
   atmosphere = Radtrans(line_species = ['H2O'])

This should produce the following output:

.. code-block:: bash
		
     Read line opacities of H2O...
    Done.

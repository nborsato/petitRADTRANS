Installation
============

Prerequisites
_____________

- A working Python installation, including the numpy package.
- A Fortran compiler (e.g. “gfortran”)


Download Code
_____________

Until the accompanying paper is published, please `contact`_ Paul Mollière in order to be added to the GitLab repository for the project. Only then will you be
able to clone it.

.. _contact: molliere@strw.leidenuniv.nl

Download petitRADTRANS from `Gitlab <https://gitlab.com/mauricemolli/petitRADTRANS.git>`_, or clone it from GitLab via

.. code-block:: bash
		
   git clone git@gitlab.com:mauricemolli/petitRADTRANS.git

.. note::
   We are working on making petitRADTRANS available via pip install and anaconda, and providing as setup script. In the meantime, please carry out the (few) steps below.

Download opacity data
_____________________

Download the `opacity and input data
<https://www.icloud.com/iclouddrive/0vEY-uMZYwPMX573_TFiR7t2A#input_data>`_
(2.3 GB), unzip them, and put the "input_data" folder into the
"petitRADTRANS" folder (i.e. the same folder where the source is, if
you clone from gitlab, this should be the
petitRADTRANS folder *in* the petitRADTRANS folder). This contains the necessary files to run petitRADTRANS, and the low resolution (:math:`\lambda/\Delta\lambda=1000`) opacity files. The high resolution (:math:`\lambda/\Delta\lambda=10^6`) opacity data is too large to be fully store online (~100 GB). Please contact us for having the high-resolution opacity files you are interested in (see :ref:`avail_opas`) made available via ftp, upon request. We are working on a more permanent solution for sharing these data.

Installation
____________

- In the terminal, enter the petitRADTRANS folder (``cd petitRADTRANS``).
- Type the following in the terminal ``chmod +x make.sh``, and press Enter.
- Type the following in the terminal ``./make.sh``, and press Enter. A lot of text will appear while the Fortran subroutines are being built. If you use Anaconda, see the first installation tip below before carrying out this step.
- Type ``ls`` in the terminal and press Enter. If everything has worked you should see three files with the “.so” extension in the folder. If you are experiencing problems see the installation tips below.
- Open the “.bash_profile” or “.bashrc” file (depending on your operating system) in your home directory. Add the following as a new line to the file (you may have to use sudo for modifying this file).

.. code-block:: bash
		
   export PYTHONPATH=Path to the folder containing petitRADTRANS/:$PYTHONPATH

.. attention::
   Don’t forget to adapt the path in the line above :) ! If you clone petitRADTRANS from gitlab, this
   should be the path of the petitRADTRANS folder *in* the petitRADTRANS folder. If you are
   uncertain what the absolute path of the folder containing the
   petitRADTRANS folder is, then switch to that folder in the
   terminal, type “pwd”, and press Enter. Don’t forget the dash “/“
   behind the path.
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

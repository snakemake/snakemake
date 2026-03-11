Setup on Windows
::::::::::::::::

If you already use Linux or MacOS X, go on with **Step 1**.

Windows Subsystem for Linux
"""""""""""""""""""""""""""

If you use Windows 10, you can set up the Windows Subsystem for Linux (`WSL`_) to natively run linux applications.
Install the WSL following the instructions in the `WSL Documentation`_. You can chose any Linux distribution available for the WSL, but the most popular and accessible one is Ubuntu.
Start the WSL and set up your account; now, you can follow the steps of our tutorial from within your Linux environment in the WSL.

Vagrant virtual machine
"""""""""""""""""""""""

If you are using a version of Windows older than 10 or if you do not wish to install the WSL, you can instead setup a Linux virtual machine (VM) with Vagrant_.
First, install Vagrant following the installation instructions in the `Vagrant Documentation`_.
Then, create a new directory you want to share with your Linux VM, for example, create a folder named ``vagrant-linux`` somewhere.
Open a command line prompt, and change into that directory.
Here, you create a 64-bit Ubuntu Linux environment with

.. code:: console

    > vagrant init hashicorp/precise64
    > vagrant up

If you decide to use a 32-bit image, you will need to download the 32-bit version of Miniconda in the next step.
The contents of the ``vagrant-linux`` folder will be shared with the virtual machine that is set up by vagrant.
You can log into the virtual machine via

.. code:: console

    > vagrant ssh

If this command tells you to install an SSH client, you can follow the instructions in this Blogpost_.
Now, you can follow the steps of our tutorial from within your Linux VM.
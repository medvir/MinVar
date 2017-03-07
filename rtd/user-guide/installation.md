## Installation (work in progress)

A condensed list of dependencies

- python 3.x
- pandas
- BLAST
- samtools
- GATK
- picard-tools
- bwa
- lofreq

### Install GATK

GATK is not (yet?) included in the instructions below because you need
to register online [here](https://software.broadinstitute.org/gatk/download/)
in order to download it. Then, install it in `/usr/local/GATK`. MinVar expects
to find the Java archive file `/usr/local/GATK/GenomeAnalysisTK.jar`.

### Setting up with Ansible

The directory [`ansible`](https://github.com/ozagordi/MinVar/tree/master/ansible)
in the project repo provides an easy option to set up MinVar with
its dependencies on a dedicated machine. The files therein define an
[ansible](http://www.ansible.com/) play for this purpose.

#### What is Ansible

Ansible is a deployment tool that allows an automatic provisioning of machines
on the cloud (it can be used on AWS, DigitalOcean, Google Cloud Platform etc.).
The user installs it on a local machine (this can be your old laptop), defines
commands in specific files and uses them to set up a remote machine.

In the following we will assume that you have installed ansible on your local
laptop and you want to set up MinVar on a remote machine running Linux Ubuntu 16.04.
The access to this machine is provided by private-public SSH key pair that must be
set up. Good instructions for this task can be found in this
[help](https://help.ubuntu.com/community/SSH/OpenSSH/Keys).

#### How to proceed

1. Install ansible on your machine (on Mac OS X you can do it via the package
   manager [brew](https://brew.sh) with `brew install ansible`),
2. clone MinVar from GitHub with `git clone https://github.com/ozagordi/MinVar.git`
   or download/unzip it,
3. move to the directory `ansible` in the cloned repository and identify the file
   `hosts`. This file contains two lines `[minvarmachine]` and a fake ip address.
   Edit the address to the one of the machine you want to setup,
4. in the same directory edit the file `setup-hosts.sh`: adapt
   `--key-file=path_to_your_private_key` to point to your personal **private** key,
5. copy your personal **public** key into `my.key`,
6. run `./setup-hosts.sh`.

#### What can go wrong

We assumed that you have an Ubuntu 16.04 available. Most of the stuff will work
on Ubuntu 14.04, but you need to switch to a newer version of Java provided by
Oracle in order to run `picard`. You can find
[here](https://www.digitalocean.com/community/tutorials/how-to-install-java-with-apt-get-on-ubuntu-16-04)
a good tutorial on this.

`UNREACHABLE!` usually means that the private/public key pair does not work.
You must be able to ssh into the remote machine with this command (edit accordingly)

    ssh -i path_to_your_private_key ubuntu@remote_ip_address

### If you don't want to use ansible

The ansible [playbook](https://github.com/ozagordi/MinVar/blob/master/ansible/setup.yml)
reads almost as plain English. You can manually copy the instructions from there
and install what you need.

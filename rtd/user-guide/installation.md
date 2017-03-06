## Installation

A condensed list of dependencies

- python 3.x
- pandas
- BLAST
- samtools
- GATK
- picard-tools
- bwa
- lofreq


### Setting up with Ansible

The directory [`ansible`](https://github.com/ozagordi/MinVar/tree/master/ansible)
in the project repo provides an easy option to set up MinVar with
its dependencies on a dedicated machine. The files therein define an
[ansible](http://www.ansible.com/) play for this purpose.

#### What is Ansible

Ansible is a deployment tool that allows an automatic provisioning of machines
on the cloud (it can be used on AWS, DigitalOcean, Google Cloud Platform etc.)
The user install ansible on their own local machine (this can be your old laptop),
defines ansible commands in specific files and uses them to set up a remote
machine.

In the following we will assume that you have installed ansible on your local
laptop and you want to set up MinVar on a remote machine running Ubuntu Linux.

https://www.digitalocean.com/community/tutorials/how-to-install-java-with-apt-get-on-ubuntu-16-04

#!/bin/bash
ansible-playbook setup.yml -i hosts --key-file=~/.ssh/oct2016.key -v # --check

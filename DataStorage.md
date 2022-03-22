# Data Storage

Goal: mount volumes ready for use on a running instance.

Three steps:

* Create volume (web interface handles this: an fdisk command)
* Format  (command line mkfs command : `sudo fdisk -l /dev/vdc`)
* Mount (two steps: mkdir data and then `sudo mount /dev/vdc`)

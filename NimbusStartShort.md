# Summary of Steps by Nimbus Instance users to achieve SSH login

Post date: 26.7.21

Author: Craig Duncan

##	Creating key pairs

###	Linux

On linux you can create key pairs with the command line: ssh-keygen

If you create these key pairs on linux, they will be named something like id_rsa (private key) and id_rsa.pub (public key).  The default location and name of these keys on a linux machine is id_rsa and they are stored in the hidden folder /.ssh/  

(To see a hidden folder you need to type ls -a instead of the usual ls at the command line)

###	Via the host

If you create the key pairs on the Pawsey machine they will have a .pem extension: e.g. yourpublickey.ppk and yourprivatekey.pem.   

The ‘pem’ file is a file that contains the private key, but has additional information to the linux RSA files.   Both types of public key files (.pem or id_rsa.pub) should be recognised and valid on the host/server machine.

###	Windows

Windows users must complete an additional step, by opening the .pem file up in putty and resaving it as a .ppk file.

See https://docs.openstack.org/horizon/latest/user/configure-access-and-security-for-instances.html#keypair-add 

## Saving the public key on the host (Pawsey server).

The public key must be saved onto the Pawsey server before logging in via SSH.  

This can be done by the instance administrator:

(a) 	uploading a public key, or 

(b)	create both public/private keys using the Pawsey computer (e.g. ‘nova’), and then saving the private key to the user’s local computer.   

##	Ensuring the local SSH program knows what the private key is and can use it

The SSH command ‘ssh add privatekeyname’ should be entered before attempting to login, so that the SSH program knows what locally saved private keys can be tested for compatibility with the public key on the host.

The private key should also have its permissions set to ‘600’ with:
chmod 600 privatekeyname

##	Logging in remotely

You need to login in to the Pawsey instance using SSH with an ‘alternative’ username to the user name you use on your own PC/laptop.   

The username for SSH depends on the operating system installed by Pawsey.  If it is a ubuntu-based instance, the chosen username will be ‘ubuntu’.   

The IP address you need to use is the floating IP address provided in the instance setup.

Then, to open the SSH connection with the default public/private key pair, you use:
ssh -Y ubuntu@xxx.xxx.xxx.xxx    (where the xxx.xxx.xxx.xxx is your IP address)

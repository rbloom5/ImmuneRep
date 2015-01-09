#!/usr/bin/python

from optparse import OptionParser
import logging
import time
import subprocess
from boto.ec2.connection import EC2Connection
import boto

def terminate_instances(ec2, key_name):
    for reservation in ec2.get_all_instances():
        for instance in reservation.instances:
            if instance.state not in ['shutting-down', 'terminated'] and instance.key_name == key_name:
                logging.info('terminating instance %s (%s)', instance.id, instance.state)
                ec2.terminate_instances([instance.id])

def print_instances(ec2):
    for reservation in ec2.get_all_instances():
        for instance in reservation.instances:
            print instance.id, instance.state, instance.ip_address

# Returns True if everything went OK 
#
def execute_command(cmd):
    p = subprocess.Popen(cmd,
                         stdin=subprocess.PIPE,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT)
    (output, error) = p.communicate()
    logging.info('%s exited; returncode: %s; stdout: %s; stderr: %s', cmd, p.returncode, repr(output), repr(error))
    return p.returncode == 0


# Create an instance and wait for it to be ssh-able
#
def run_instance(ec2, image_id, instance_type, key_name, key_file, timeout=300):
	start_time = time.time()
	reservation = ec2.run_instances(image_id=image_id, instance_type=instance_type, key_name=key_name)
	instance = reservation.instances[0]
	logging.info('started instance %s', instance.id)
	while time.time() - start_time < timeout:
	    if instance.update() == 'running':
	        instance = ec2.get_all_instances([instance.id])[0].instances[0]
	        logging.info('instance %s now running at %s', instance.id, instance.ip_address)
	        while time.time() - start_time < timeout:
	            cmd = ['ssh', '-o', 'StrictHostKeyChecking=no', '-i', key_file, 'ec2-user@' + instance.ip_address, 'echo']
	            if execute_command(cmd):
	                logging.info('instance %s responding %s', instance.id, instance.ip_address)
	                return instance.ip_address
	            time.sleep(2)
	        pass
	    time.sleep(2)

def main():
    usage = "usage: %prog [options] remote-command [remote-args...]"
    parser = OptionParser(usage=usage)

    parser.add_option('--aws-access-key-id', help='AWS access key id (required)')
    parser.add_option('--aws-secret-access-key', help='AWS secret access key (required)')
    parser.add_option('--instance-type', help='The type of instance, e.g. m1.small (required); see http://aws.amazon.com/ec2/instance-types/')
    parser.add_option('--ami-id', help='The id of the machine image (required); use ami-84db39ed for 32-bit and ami-86db39ef for 64-bit Fedora Linux 8')
    parser.add_option('--key-name', help='Name of the key pair (required)')
    parser.add_option('--key-file', help='Filepath of the key pair (required)')
    parser.add_option('-i', '--input-file', action='append', help='A file to trasfer to the instance before running the remote command')
    parser.add_option('-o', '--output-file', action='append', help='A file to trasfer from the instance after running the remote command')
    parser.add_option('--log', help='A file to append progress / error information to')
    
    (options, args) = parser.parse_args()
    print options.ami_id
    print args
    
    required_options = ['aws_access_key_id', 'aws_secret_access_key', 'instance_type', 'ami_id', 'key_name', 'key_file']
    for option in required_options:
        if not options.__dict__[option]:
            print 'missing required option', option
            print 
            parser.print_help()
            exit(1)
    if not args:
        print 'missing remote-command'
        print 
        parser.print_help()
        exit(1)

    if options.log:
        logging.basicConfig(level=logging.INFO,
                            format='%(asctime)s %(levelname)s %(message)s',
                            filename=options.log,
                            filemode='a')
        logging.getLogger('boto').setLevel(level=logging.ERROR)

    ec2 = boto.ec2.connect_to_region('us-west-1', aws_access_key_id=options.aws_access_key_id, aws_secret_access_key=options.aws_secret_access_key)
    terminate_instances(ec2, options.key_name)
    try:
        instance_ip = run_instance(ec2, options.ami_id, options.instance_type, options.key_name, options.key_file)
        if instance_ip:
            remote_user_host = 'ec2-user@' + instance_ip
            if options.input_file:
                for f in options.input_file:
                    cmd = ['scp', '-o', 'StrictHostKeyChecking=no', '-i', options.key_file, f, remote_user_host + ':./']
                    returncode = execute_command(cmd)                
            cmd = ['ssh', '-o', 'StrictHostKeyChecking=no', '-i', options.key_file, remote_user_host]
            cmd.extend(args)
            returncode = execute_command(cmd)
        else:
            raise Error('unable to connect to instance')
    finally:
        terminate_instances(ec2, options.key_name)
        
main()


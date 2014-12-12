#!/usr/bin/env python

import boto
import boto.ec2
from boto.s3.connection import S3Connection
from boto.s3.connection import Location
from boto.s3.key import Key


def S3_connect():
	pass


def start_EC2(instanceID = 'i-4f290d85'):
	#will start up an EC2 instance.  Default ID is the 'vdj-pipeline-server'
	conn = boto.ec2.connect_to_region('us-west-1', \
		aws_access_key_id = 'AKIAJ2TEUHQV2LHU7XQQ', \
		aws_secret_access_key = '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')

	conn.start_instances('i-4f290d85')


def get_instance_ids():
	conn = boto.ec2.connect_to_region('us-west-1', \
	aws_access_key_id = 'AKIAJ2TEUHQV2LHU7XQQ', \
	aws_secret_access_key = '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')
	instances = []
	for r in conn.get_all_reservations():
		instances.append(r.instances[0])

	return instances
import json
from boto.s3.connection import S3Connection
from boto.s3.key import Key

def update_default_groups(name, files):
    conn = S3Connection('AKIAJ2TEUHQV2LHU7XQQ','VKzoYINBlvZi5uiAIeAnJG5fgLedQPFrmMCpSfBp')
    info_buck = conn.get_bucket('info-docs')
    info_buck.get_key("ExperimentClassDefaultGroups.txt").get_contents_to_filename("ExperimentClassDefaultGroups.txt")
    with open("ExperimentClassDefaultGroups.txt") as default_groups_location:
        defaultgroups=json.load(default_groups_location)
        
    defaultgroups[name] = files
    
    with open("ExperimentClassDefaultGroups.txt", 'w') as f:
        json.dump(defaultgroups, f)
        
    info_buck.get_key("ExperimentClassDefaultGroups.txt").set_contents_to_filename("ExperimentClassDefaultGroups.txt")

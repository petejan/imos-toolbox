#!/usr/bin/env python3
""" Build the imos stand-alone toolbox GUI using only the Matlab compiler.

This script ignores untracked and testfiles within the repository and weakly
attach a string version to the binary filename.
The string is related to the repo state.

Usage:
  build.py  --arch=<architecture> [--root_path=<imostoolboxpath> --mcc_path=<mccpath> --dist_path=<distpath>]

Options:
  -h --help     Show this screen.
  --version     Show version.
  --mcc_path=<mccpath>           The Matlab runtime Path
  --root_path=<imostoolboxpath>  The repository root path
  --dist_path=<distpath>  Where Compiled items should be stored
  --arch=<architecture>   One of win64,glnxa64,maci64
"""

import os
import subprocess as sp
from getpass import getuser
from pathlib import Path
from shutil import copy
from typing import List

from docopt import docopt
from git import Repo

VALID_ARCHS = ['glnxa64', 'win64', 'maci64']
ARCH_NAMES = {
    'glnxa64': 'Linux64',
    'maci64': 'Mac64',
    'win64': 'Win64',
}
VERSION_FILE = '.standalone_canonical_version'


def run(x: str):
    proc = sp.run(x, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
    return proc.stdout.decode('utf-8').strip(), proc.stderr.decode(
        'utf-8').strip()


def create_java_call_sig_compile(root_path: str) -> str:
    java_path = os.path.join(root_path, 'Java')
    return f"cd {java_path} && ant install && cd {root_path}"


def git_info(root_path: str) -> dict:
    """ Return git information
    info['username']
    info['is_dirty']
    info['branch']
    info['tag']
    info['modified_files']
    info['staged_files']
    info['untracked_files']
    info['is_tag_release']
    info['is_master']
    info['is_clean']
    info['is_official_release']
    info['version']
    """
    repo = Repo(root_path)
    if repo.bare:
        raise TypeError(f"{root_path} is not a git repository")

    info = {}
    info['username'] = getuser()
    info['is_dirty'] = repo.is_dirty()
    info['branch'] = str(repo.active_branch)
    info['tag'] = str(repo.git.describe('--tags'))
    info['modified_files'] = [x.a_path for x in repo.index.diff(None)]
    info['staged_files'] = [x.a_path for x in repo.index.diff("HEAD")]
    info['untracked_files'] = repo.untracked_files

    info['is_tag_release'] = len(info['tag'].split('-')) < 3
    info['is_master'] = info['branch'] == 'master'
    info['is_clean'] = not info['is_dirty'] and not info['staged_files']
    info['is_official_release'] = info['is_tag_release'] and info['is_clean']

    version = ''
    if info['is_official_release']:
        version += info['tag']
    else:
        version += info['branch']
        version += '-' + info['username']
        version += '-' + info['tag']
        if info['is_dirty']:
            version += '-dirty'

    info['version'] = version
    return info


def find_files(folder: str, ftype: str) -> list:
    """ Find all the files within the repository that are not testfiles """
    for file in Path(folder).glob(os.path.join('**/', ftype)):
        filepath = file.as_posix()
        if 'data/testfiles' not in filepath and 'tmpbuild.m' not in filepath:
            yield filepath


def get_required_files(root_path: str, info: dict) -> (list, list):
    """ obtain the mfiles and matfiles for binary compilation """
    allmfiles = list(find_files(root_path, '*.m'))
    allmatfiles = list(find_files(root_path, '*.mat'))
    if info['is_dirty']:
        print(f"Warning: {root_path} repository is dirty")
    if info['modified_files']:
        print(f"Warning: {root_path} got modified files not staged")
    if info['staged_files']:
        print(f"Warning: {root_path} got staged files not commited")
    if info['untracked_files']:
        print(
            f"Warning: {root_path} got untracked files - Ignoring them all...")
        for ufile in info['untracked_files']:
            if ufile in allmfiles:
                allmfiles.remove(ufile)
            if ufile in allmatfiles:
                allmatfiles.remove(ufile)

    mind = [
        ind for ind, file in enumerate(allmfiles) if 'imosToolbox.m' in file
    ]
    mentry = allmfiles.pop(mind[0])
    allmfiles = [mentry] + allmfiles
    return allmfiles, allmatfiles


def create_mcc_call_sig(mcc_path: str, root_path: str, dist_path: str,
                        mfiles: List[str], matfiles: List[str],
                        output_name: str, arch: str) -> List[str]:
    standalone_flags = '-m'
    verbose_flag = '-v'
    outdir_flag = f"-d {dist_path}"
    clearpath_flag = '-N'
    outname_flag = f"-o {output_name}"
    verbose_flag = '-v'
    warning_flag = '-w enable'
    mfiles_str = ' '.join([f"{file}" for file in mfiles])
    matfiles_str = ' '.join([f"-a {file}" for file in matfiles])

    mcc_arguments = ' '.join([
        standalone_flags, verbose_flag, outdir_flag, clearpath_flag,
        outname_flag, verbose_flag, warning_flag, mfiles_str, matfiles_str
    ])

    if arch == 'win64':
        tmpscript = os.path.join(root_path, 'tmpbuild.m')
        with open(tmpscript,
                  'w') as tmpfile:  # overcome cmd.exe limitation of characters
            tmpfile.write(f"eval('{mcc_path} {mcc_arguments}')")
        external_call = f"matlab -nodesktop -nosplash -wait -r \"run('{tmpscript}');exit\""
        rm_call = f"del {tmpscript}"
        return [external_call, rm_call]
    else:
        return [
            ' '.join([
                mcc_path, standalone_flags, verbose_flag, outdir_flag,
                clearpath_flag, outname_flag, verbose_flag, warning_flag,
                mfiles_str, matfiles_str
            ])
        ]


def get_args(args: dict) -> (str, str, str, str):
    """process the cmdline arguments to return the root_path, dist_path, mcc binary, and architecture"""
    if not args['--root_path']:
        root_path = os.getcwd()
    else:
        root_path = args['--root_path']

    if not args['--dist_path']:
        dist_path = os.path.join(root_path, 'dist')
        exists = os.path.lexists(dist_path)
        if not exists:
            os.mkdir(dist_path)
    else:
        dist_path = args['--dist_path']

    if not args['--mcc_path']:
        mcc_path = 'mcc'
    else:
        mcc_path = args['--mcc_path']

    if args['--arch'] in VALID_ARCHS:
        ind = VALID_ARCHS.index(args['--arch'])
        arch = VALID_ARCHS[ind]
    else:
        raise Exception(
            f"{args['--arch']} is not one of the valid architectures {VALID_ARCHS}"
        )
    return root_path, dist_path, mcc_path, arch


def write_version(afile: str, version: str) -> None:
    with open(afile, 'w') as file:
        file.writelines([version + '\n'])


if __name__ == '__main__':
    args = docopt(__doc__, version='0.1')

    root_path, dist_path, mcc_path, arch = get_args(args)
    repo_info = git_info(root_path)

    print("Starting build process.")
    print(f"Marking {repo_info['version']} as the standalone version")
    write_version(VERSION_FILE, repo_info['version'])

    # output name of binary is restricted in mcc
    tmp_name = f"imosToolbox_{ARCH_NAMES[arch]}"

    print("Gathering files to be included")
    mfiles, matfiles = get_required_files(root_path, repo_info)
    java_call = create_java_call_sig_compile(root_path)

    mcc_call = create_mcc_call_sig(mcc_path, root_path, dist_path, mfiles,
                                   matfiles, tmp_name, arch)

    print("Current repo information:")
    print(repo_info)

    print(f"Calling {java_call}..")
    stdout, stderr = run(java_call)
    if stderr:
        __import__('pdb').set_trace()
        raise Exception(f"{jcall} failed")

    print(f"Calling {mcc_call}...")
    for mcall in mcc_call:
        stdout, stderr = run(mcall)
        if stderr:
            raise Exception(f"{mcall} failed")

    print(f"The toolbox architecture is {ARCH_NAMES[arch]}.")
    print(f"Git version information string is {repo_info['version']}")

    print(f"Updating binary at {root_path}....")

    # mcc append .exe to end of file
    if arch == 'win64':
        copy(os.path.join(dist_path, tmp_name + '.exe'),
             os.path.join(root_path, tmp_name + '.exe'))
    else:
        copy(os.path.join(dist_path, tmp_name),
             os.path.join(root_path, tmp_name + '.bin'))

    print("Build Finished.")

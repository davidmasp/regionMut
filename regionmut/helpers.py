import subprocess
import re
import os
import tempfile
import sys


def find_r_installed():
    """
    Find R installed.
    """
    outpt = subprocess.run(['which', 'Rscript'], check=True, capture_output=True)
    path_to_r = outpt.stdout.decode("utf-8").strip()
    return path_to_r

def find_r_version(path_r):
    """
    Find which version has R installed
    """
    cmd = [path_r, '--version']
    outpt = subprocess.run(cmd, check=True, capture_output=True)
    version_out = outpt.stderr.decode("utf-8").strip()
    vout = re.findall(r'(\d+\.\d+\.\d+)', version_out)[0]
    return vout

def find_r_lib(path_r):
    """
    pass
    """

def get_exec_path(script, pkg):
    """
    Get path to executable
    """
    pathR = find_r_installed()
    r_shibang = "#!{}".format(pathR)
    r_code = 'system.file("{}", package = "{}")'.format(script, pkg)
    tmp_script = "\n".join([
        r_shibang,
        "",
        r_code
    ])
    tmp_file = tempfile.NamedTemporaryFile(mode = "w", delete = False, dir = ".")
    tmp_file.write(tmp_script)
    tmp_file.close()
    scrp_out = subprocess.run([pathR, tmp_file.name], check=True, capture_output=True)
    os.remove(tmp_file.name)
    tout = scrp_out.stdout.decode("utf-8").strip()
    # from https://stackoverflow.com/a/9085630
    pattern = r'"([A-Za-z0-9_\./\\-]*)"'
    script_path = re.findall(pattern, tout)[0]
    return script_path

def get_exec_path_mode(mode, pkg):
    ### this needs to be standardized
    # exec/pkg_mode.R
    script = "exec/{}_{}.R".format(pkg, mode)
    return get_exec_path(script, pkg)


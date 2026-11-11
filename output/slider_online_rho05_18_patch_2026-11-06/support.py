"""Atomic code-installation helpers; no fitting or scheduling."""

import hashlib, importlib.util, json, os, shutil, sys

from pathlib import Path

sys.dont_write_bytecode=True

HERE=Path(__file__).resolve().parent

def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''): h.update(chunk)
    return h.hexdigest()

def read(path): return json.loads(Path(path).read_text())

def encoded(obj): return (json.dumps(obj,indent=2,sort_keys=True)+'\n').encode()

def write(path,data):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_name(path.name+'.repair_tmp')
    with temp.open('wb') as f:f.write(data);f.flush();os.fsync(f.fileno())
    os.replace(str(temp),str(path))

def atomic(path,obj):write(path,encoded(obj))

def inside(path,root):
    path=Path(path)
    if path.is_symlink() or root.resolve() not in path.resolve().parents:
        raise RuntimeError('Unsafe repair path: '+str(path))
    return path

def move(source,destination,root):
    inside(source,root);inside(destination,root)
    if destination.exists():raise RuntimeError('Archive already exists: '+str(destination))
    destination.parent.mkdir(parents=True,exist_ok=True)
    os.rename(str(source),str(destination))

def verify_manifest(code,manifest):
    for line in manifest.decode().splitlines():
        expected,name=line.split('  ',1)
        if digest(inside(code/name,code))!=expected:raise RuntimeError('Code changed: '+name)

def runtime():
    spec=importlib.util.spec_from_file_location('submission_repair_runtime',HERE/'runtime/minibatch.py')
    mb=importlib.util.module_from_spec(spec);sys.modules[spec.name]=mb;spec.loader.exec_module(mb);return mb

def patch_code(branch,state,old,new):
    code=branch/'code';backup=state/'code.before'
    if code.exists():
        current=(code/'SHA256SUMS').read_bytes()
        if current==new:verify_manifest(code,new);return
        if current!=old:raise RuntimeError('Unexpected installed code')
        verify_manifest(code,old)
    elif not backup.exists():raise RuntimeError('Installed code and backup both missing')
    stage=state/'code.next'
    if not stage.exists():shutil.copytree(str(HERE/'runtime'),str(stage),ignore=shutil.ignore_patterns('__pycache__'))
    verify_manifest(stage,new)
    if code.exists():move(code,backup,branch)
    move(stage,code,branch)

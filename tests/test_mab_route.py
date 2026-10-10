import os
import random
import shutil
import subprocess
import sys
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[1]

def test_mab_route_calls_snp_and_preserves_read_positions(tmp_path):
    if any(shutil.which(t) is None for t in ('bwa', 'samtools', 'java', 'perl')):
        pytest.skip('alignment and calling tools are required')
    rng=random.Random(13)
    reference=''.join(rng.choice('ACGT') for _ in range(1500))
    alternative=next(b for b in 'ACGT' if b != reference[499])
    strain=reference[:499]+alternative+reference[500:]
    ref=tmp_path/'ref.fa';ref.write_text('>chr\n'+reference+'\n')
    def fastq(path, sequences):
        path.write_text(''.join(f'@r{i}\n{s}\n+\n'+ 'I'*len(s)+'\n' for i,s in enumerate(sequences)))
    r1=tmp_path/'r1.fq';r2=tmp_path/'r2.fq';single=tmp_path/'single.fq'
    complement=str.maketrans('ACGT','TGCA')
    fastq(r1,[strain[400:550]]*12)
    fastq(r2,[strain[600:750].translate(complement)[::-1]]*12)
    fastq(single,[strain[425:575]]*4)
    bam=tmp_path/'test.bam';helper=ROOT/'src/scripts/NTM_bwamem.sh'
    subprocess.run(['bash',str(helper),'align',str(ref),str(ref),str(r1),str(r2),str(single),str(bam),'test','1','1'],check=True,cwd=ROOT,capture_output=True)
    subprocess.run(['samtools','faidx',str(ref)],check=True)
    mask=tmp_path/'mask';mask.write_text('')
    out=tmp_path/'variants'
    subprocess.run(['bash',str(helper),'call',str(ref),str(bam),str(ROOT/'src/scripts/VarScan.v2.3.9.jar'),str(mask),str(out),'test',sys.executable,'src.cli'],check=True,cwd=ROOT,capture_output=True)
    assert f'500\t{reference[499]}\t{alternative}' in (out/'test.snp').read_text()
    records=(out/'test.forup').read_text().splitlines()
    assert records
    position_field=records[0].split('\t')[-2]
    assert all(n.isdigit() for n in position_field.split(','))
    assert (out/'test.cns').stat().st_size > 0
    assert not list(out.glob('.mab-call.*'))

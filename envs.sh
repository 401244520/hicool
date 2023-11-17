conda create -n hicool python=3.10
conda activate hicool
pip install -i https://pypi.tuna.tsinghua.edu.cn/simple cooltools
pip install pyBigWig tqdm
conda install dash 
pip install -e .
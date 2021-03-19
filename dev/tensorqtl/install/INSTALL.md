### Setup CUDA drivers and PyTorch on GCP

Launch a new instance configured with Ubuntu 20.04 LTS and a GPU, clone this repository, and run the following:
#### Install CUDA
```bash
sudo ./install_cuda.sh
sudo reboot
# verify
nvidia-smi
```

#### Install Python 3
```bash
sudo apt update
sudo apt install -y python3-pip python3-dev python3-virtualenv
pip3 install --upgrade pip virtualenv setuptools
virtualenv venv
source venv/bin/activate
pip3 install -r requirements.txt
# verify
python -c "import torch; print(torch.__version__); print('CUDA available: {} ({})'.format(torch.cuda.is_available(), torch.cuda.get_device_name(torch.cuda.current_device())))"
```

#### Install rmate (optional)
```bash
sudo apt install -y ruby
mkdir ~/bin
curl -Lo ~/bin/rmate https://raw.githubusercontent.com/textmate/rmate/master/bin/rmate
chmod a+x ~/bin/rmate
echo 'export RMATE_PORT=${rmate_port}' >> ~/.bashrc
```

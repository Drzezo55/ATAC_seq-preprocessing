
# Download Picard JAR
mkdir -p ~/tools/picard
cd ~/tools/picard
wget https://github.com/broadinstitute/picard/releases/download/3.1.1/picard.jar
# Optional: Create a shell alias in your ~/.bashrc or ~/.zshrc
echo "alias picard='java -jar ~/tools/picard/picard.jar'" >> ~/.bashrc
source ~/.bashrc
# Clone and install HOMER
cd ~
git clone https://github.com/timydaley/homer.git
cd homer
perl configureHomer.pl -install
# Add HOMER to PATH
echo 'export PATH="$HOME/homer/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc

# Install mm10 genome
perl configureHomer.pl -install mm10
# install other tools via conda

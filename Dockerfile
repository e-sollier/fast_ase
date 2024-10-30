FROM rust:1.82.0
RUN apt-get update && apt-get install -y cmake libncurses5-dev libbz2-dev liblzma-dev clang
RUN wget https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2 &&  tar -xjf htslib-1.21.tar.bz2 && \
 rm htslib-1.21.tar.bz2 && cd htslib-1.21 && ./configure && make && cp bgzip /bin && cp tabix /bin && cp htsfile /bin && cd ..
COPY . .
RUN cargo install --path .


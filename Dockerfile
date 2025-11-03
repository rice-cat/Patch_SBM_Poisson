# Use official deal.II Docker image
FROM dealii/dealii:master-jammy

# Set working directory
WORKDIR /workspace

# Copy the project files
COPY . /workspace/

# Create build directory and build the project
RUN mkdir -p build && \
    cd build && \
    cmake .. && \
    make -j$(nproc)

# Set the default command to run tests
CMD ["bash"]

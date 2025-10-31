# Use a Python base image
FROM python:3.10.12-slim

# Set working directory
WORKDIR /app

# Copy only requirements first, install dependencies
COPY requirements.txt .

# Install dependencies & clean up
RUN apt update && \
    apt install -y bcftools && \
    pip install --no-cache-dir -r requirements.txt && \
    rm -rf /var/lib/apt/lists/* /root/.cache/pip

# Copy application files
COPY app/ /app/

# Set appropriate permissions
RUN chmod -R 777 /app

# Expose the port for Streamlit
EXPOSE 8501

# Create .streamlit directory and copy config/credentials files
RUN mkdir -p ~/.streamlit
COPY app/config.toml ~/.streamlit/config.toml
COPY app/credentials.toml ~/.streamlit/credentials.toml

# Set the default command
ENTRYPOINT ["streamlit", "run"]
CMD ["application.py"]

docker run -d --name clickhouse-server --ulimit nofile=262144:262144 -p 19000:9000 -p 18123:8123 -p 19009:9009 yandex/clickhouse-server

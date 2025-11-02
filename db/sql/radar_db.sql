
CREATE TABLE <table_name>
(
    id UInt64,
    code UInt64,
    lon Float32,
    lat Float32,
    alt Float32,
    density Float32,
    timestamp Float32  -- do not use time stamp since clickhouse does not support millisecond
)
ENGINE = MergeTree
ORDER BY (id, timestamp, code)
SETTINGS index_granularity = 8192

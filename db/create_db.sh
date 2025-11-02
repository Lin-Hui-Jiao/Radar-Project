db_endpoint="${1:-"http://localhost:18123"}"
db_table="${2:-"radar"}"

echo "DROP TABLE $2" | curl "$db_endpoint" --data-binary @-
cat db/sql/radar_db.sql | sed "s/<table_name>/$2/" | curl "$db_endpoint" --data-binary @-

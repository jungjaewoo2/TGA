#define SQL_STMT_DELETE_ALL \
    "PRAGMA writable_schema = 1;" \
    "delete from sqlite_master where type in ('table', 'index', 'trigger');" \
    "PRAGMA writable_schema = 0;" \
    "VACUUM;"

#define SQL_STMT_DROP_TABLE_IF_EXISTS \
    "DROP TABLE IF EXISTS ?;"

int sql_delete_all_tables(Rfile * rf);
int sql_drop_table_if_exists(Rfile * rf, char *table);

#include "rlab.h"
#include <sqlite3.h>
#include <stdio.h>

#include "r_sqlite3.h"

#undef  THIS_FUNCTION
#define THIS_FUNCTION "sql_delete_all_tables"
int sql_delete_all_tables(Rfile * rf)
{
  sqlite3_stmt *res;

  if (!rf->sql_db)
    return 1;

  int rc = sqlite3_prepare_v2(rf->sql_db, SQL_STMT_DELETE_ALL, -1, &res, 0);
  if (rc != SQLITE_OK)
    return 1;

  rc = sqlite3_step(res);
  if (rc != SQLITE_DONE)
  {
    printf(THIS_FUNCTION ": 'sqlite3_step' returned code %i, %s", rc,
           sqlite3_errmsg(rf->sql_db));
  }
  sqlite3_finalize(res);
  return 0;
}
#undef  THIS_FUNCTION
#define THIS_FUNCTION "sql_drop_table_if_exists"
int sql_drop_table_if_exists(Rfile * rf, char *table)
{
  int rval=0;
  sqlite3_stmt *res;

  if (!rf->sql_db)
    return 1;

  int rc = sqlite3_prepare_v2(rf->sql_db, SQL_STMT_DROP_TABLE_IF_EXISTS, -1, &res, 0);
  if (rc != SQLITE_OK)
  {
    rval=1;
    goto _exit_func;
  }

  rc = sqlite3_bind_text(res, 1, table, -1, SQLITE_TRANSIENT);
  if( rc != SQLITE_OK) {
    printf(THIS_FUNCTION ": 'sqlite3_bind_text' returned error %s (%i)!\n",
           sqlite3_errmsg(rf->sql_db), rc);
    rval=1;
    goto _exit_func;
  }

  rc = sqlite3_step(res);
  if (rc != SQLITE_DONE)
  {
    printf(THIS_FUNCTION ": 'sqlite3_step' returned error %s (%i)!\n",
           sqlite3_errmsg(rf->sql_db), rc);
    rval=1;
    goto _exit_func;
  }

  sqlite3_finalize(res);

_exit_func:
  return rval;
}

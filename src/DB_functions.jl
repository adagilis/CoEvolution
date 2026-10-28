using DuckDB
using TidierDB

function createERCTable(db,table) 
    con = DBInterface.connect(DuckDB.DB(db))
    DBInterface.execute(con, "CREATE TABLE $table (i VARCHAR, j VARCHAR,n_edges DOUBLE,r DOUBLE,fERC DOUBLE,pval DOUBLE,shared_tips INT,gid VARCHAR,chr_combo VARCHAR,dist BIGINT,same_chr BOOL)")
    
end
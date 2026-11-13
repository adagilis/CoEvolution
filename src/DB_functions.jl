using DuckDB
using TidierDB

function createERCTable(db,table) 
    con = DBInterface.connect(DuckDB.DB(db))
    DBInterface.execute(con, "CREATE TABLE $table (i VARCHAR, j VARCHAR,n_edges DOUBLE,r DOUBLE,fERC DOUBLE,pval DOUBLE,shared_tips INT,gid VARCHAR,chr_combo VARCHAR,dist BIGINT,same_chr BOOL)")
    DBInterface.close!(con)
end

function spec_table(db,table)
    con = DBInterface.connect(DuckDB.DB(db))
    return(db_table(con,table))
end

function mean_gene(db,table,gene)
    tab = spec_table(db,table)
    vals = @eval @chain $tab begin;
        @filter(i == $(String(gene)) || j == $(String(gene)) )
        @select(fERC)
        @collect
    end
    return(mean(vals[:,1]))
end

function all_mean_genes(db, table)
    query = """
    SELECT gene, AVG(fERC) AS mean_fERC
    FROM (
        SELECT i AS gene, fERC FROM $table
        UNION ALL
        SELECT j AS gene, fERC FROM $table
    )
    GROUP BY gene
    """
    return DataFrame(DBInterface.execute(db, query))
end
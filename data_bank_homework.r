library("RODBC")


# homework 2

pro_str <- data.frame(
    STR_ID = 1:5,
    NAME = c("first pro1", "second pro1", "third pro1", "forth pro1", "fifth pro1"),
    PDB_ACCESSION = c("123bc3", "rg345", "fgh345252", "fr43552", "55463ggf"),
    RESOLUTION = c(234, 56, 786, 876, 98)
)

print(pro_str)
print(pro_str$STR_ID)
print(pro_str$PDB_ACCESSION)


pro_seq <- data.frame(
    SEQ_ID = 11:15,
    NAME = c("second pro2", "forth pro2", "first pro2", "fifth pro2", "third pro2"),
    UNIPORT_ACCESSION = c("un23", "un65", "un76", "un00", "un12"),
    STR_ID = c(2, 4, 1, 5, 3)
)

print(pro_seq)
print(pro_seq$STR_ID)


# homework 4

connection <- odbcConnect(dsn="mixtures", believeNRows=FALSE)
data <- sqlQuery(connection, "SELECT * FROM MIXTURE")

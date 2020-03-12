#'@importFrom Biobase package.version
#'
format_argg <- function(x){
txdb_info <- x$txdb
txdb_info[1] <- ""
header_txdb <- c("########################################################",
                 "#              Transcript Annotation Info:             #",
                 "########################################################")
txdb_info <- c(header_txdb, txdb_info)
bsgenome_info <- x$bsgenome
bsgenome_info[1] <- ""
if(length(bsgenome_info) == 1) bsgenome_info = "Missing Genome Sequence Infomation"
header_bsgenome <- c("########################################################",
                     "#                Genome Reference Info:                #",
                     "########################################################")
bsgenome_info <- c(header_bsgenome, bsgenome_info)
bsgenome_info <- grep("seqnames()|given",bsgenome_info,invert = TRUE,value = TRUE)
package_version_info <- paste0("exomePeak2 Version: ", package.version("exomePeak2") )
x <- x[!names(x) %in% c("txdb","bsgenome")]
argument_info <- paste0(names(x), " = " ,x)
argument_info <- gsub("\\[1\\] ","",argument_info)
header_arguments <- c("########################################################",
                      "#                 exomePeak2 Run Info                  #",
                      "########################################################")
argument_info <- c(header_arguments, "", argument_info)

return(c(argument_info,"",txdb_info,"",bsgenome_info,"",
         "--------------------------------------------------------",package_version_info))
}

########################################################
#                 exomePeak2 Run Info                  #
########################################################


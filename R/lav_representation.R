# user visible function to add 'matrix' entries in the parameter table
lavMatrixRepresentation <- function(partable, representation = "LISREL",
                                    add.attributes = FALSE,
                                    as.data.frame. = TRUE) {

    # check parameter table
    partable <- lav_partable_complete(partable)

    # get model matrices
    if(representation == "LISREL") {
        REP <- representation.LISREL(partable, target = NULL,
                                     extra = add.attributes)
    } else {
        stop("psindex ERROR: only representation \"LISREL\" has been implemented.")
    }

    partable$mat <- REP$mat
    partable$row <- REP$row
    partable$col <- REP$col

    if(as.data.frame.) {
        partable <- as.data.frame(partable, stringsAsFactors=FALSE)
        class(partable) <- c("psindex.data.frame", "data.frame")
    }

    if(add.attributes) {
        attr(partable, "ov.dummy.names.nox") <- attr(REP, "ov.dummy.names.nox")
        attr(partable, "ov.dummy.names.x")   <- attr(REP, "ov.dummy.names.x")
        attr(partable, "mmNames")  <- attr(REP, "mmNames")
        attr(partable, "mmNumber") <- attr(REP, "mmNumber")
        attr(partable, "mmRows")   <- attr(REP, "mmRows")
        attr(partable, "mmCols")   <- attr(REP, "mmCols")
        attr(partable, "mmDimNames")  <- attr(REP, "mmDimNames")
        attr(partable, "mmSymmetric") <- attr(REP, "mmSymmetric")
    }

    partable
}


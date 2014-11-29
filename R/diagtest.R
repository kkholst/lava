##' @export
diagtest <- function(table,positive=1,...) {
    if (!is.table(table) && (is.matrix(table) || is.data.frame(table))) {
        table <- base::table(table[,1],table[,2])
        names(dimnames(table)) <- colnames(table)[1:2]
    }
    if (!is.table(table) || nrow(table)!=2 || ncol(table)!=2) stop("2x2 table expected")
    if (positive==2) {
        nam <- names(dimnames(table))
        table <- as.table(table[2:1,2:1,drop=FALSE])
        names(dimnames(table)) <- nam
    }
    M <- multinomial(table)
    diag <- function(p,...) {
        list(Prevalence=(p[1]+p[2]),
             Sensitivity=p[1]/(p[1]+p[2]), # Prob test + | given (true) disease (True positive rate)
             Specificity=p[4]/(1-p[1]-p[2]), # Prob test - | given no disease (True negative rate)
             PositivePredictiveValue=p[1]/(p[1]+p[3]), # Prob disease | test +
             NegativePredictiveValue=p[4]/(1-p[1]-p[3]) # Prob disease free | test .
             )}    
    list(table=table,multinomial=M,estimate=estimate(M,diag))
}

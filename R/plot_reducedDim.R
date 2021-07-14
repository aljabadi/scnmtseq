#' Plot reduced dimensions for SCE object
#'
#' @param sce A SCE object
#' @param reducedDim Name of the reduced dimension from \code{names(sce@int_colData$reducedDims)}
#' @param comps Integer vector of length 2. The components to use.
#' @param colBy Character or integer. colData to use to colour the cells. From \code{names(colData(sce))}
#' @param facet_by Character or integer. colData to use to facet the plots. From \code{names(colData(sce))}
#' @param cols Named character of colours to use. If NULL default values are generated.
#' @param axis.title Title prefix of the 2 axes, e.g. 'PC' or 'UMAP' etc.
#' @param ... Additional parameters passed to \code{\link[ggplot2]{geom_point}}
#'
#' @return A ggplot2 object
#' @export
#'
plot_reducedDim <-
    function(sce,
             reducedDim = 'PCA',
             comps = c(1, 2),
             colBy = 'Batch',
             facet_by = NULL,
             cols = NULL,
             axis.title = 'PC',
             ...)
    {
        df <- sce@int_colData$reducedDims[[reducedDim]]
        df <- df[, comps]
        df <- as.data.frame(df)
        colnames(df) <- paste0(axis.title, comps)
        pcs <- colnames(df)
        df$group <- colData(sce)[, colBy]
        df$facet <- colData(sce)[, facet_by]

        if (is.null(cols))
        {
            n_cols <- length(unique(df$group))
            hues = seq(15, 375, length = n_cols + 1)
            cols <- grDevices::hcl(h = hues, l = 65, c = 100)[seq_len(n_cols)]
        }

        p <-
            ggplot(df, aes_string(pcs[1], pcs[2])) + geom_point(aes(col = group), ...) +
            scale_colour_manual(values = cols) +
            theme_classic() + labs(col = colBy)
        if (!is.null(facet_by))
        {
            p <- p + facet_wrap(. ~ facet, scales = 'free')
        }
        p
    }

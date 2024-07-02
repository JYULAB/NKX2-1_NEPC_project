library(cowplot)
library(patchwork)
# Meta highligh plot ----
# added order to list of params
#changed DimPlot_scCustom_man to DimPlot_scCustom_man_man, and cahnged order from true to order
Meta_Highlight_Plot_man <- function(
  seurat_object,
  meta_data_column,
  meta_data_highlight,
  highlight_color = NULL,
  background_color = "lightgray",
  pt.size = NULL,
  aspect_ratio = NULL,
  figure_plot = FALSE,
  raster = NULL,
  raster.dpi = c(512, 512),
  label = FALSE,
  split.by = NULL,
  split_seurat = FALSE,
  ggplot_default_colors = FALSE,
  order = NULL,
  ...
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)
  
  # Check meta data
  good_meta_data_column <- Meta_Present(object = seurat_object, meta_col_names = meta_data_column, omit_warn = FALSE, print_msg = FALSE, return_none = TRUE)[[1]]
  
  # stop if none found
  if (length(x = good_meta_data_column) == 0) {
    cli_abort(message = c("{.code meta_data_column} was not found.",
                          "i" = "No column found in object meta.data named: {.val {meta_data_column}}.")
    )
  }
  
  # Check that meta data is factor or character
  accepted_meta_types <- c("factor", "character", "logical")
  
  if (!class(x = seurat_object@meta.data[[good_meta_data_column]]) %in% accepted_meta_types) {
    cli_abort(message = c("The {.code good_meta_data_column}: {.field {good_meta_data_column}} is of class: {.val {class(x = seurat_object@meta.data[[good_meta_data_column]])}}.",
                          "i" = "Meta data variables must be of classes: factor, character, or logical to be used with {.code Meta_Highlight_Plot()}.")
    )
  }
  
  # Check meta_data_highlight
  meta_var_list <- as.character(x = unique(x = seurat_object@meta.data[, good_meta_data_column]))
  
  # Check good and bad highlight values
  bad_meta_highlight <- meta_var_list[!meta_var_list %in% meta_data_highlight]
  found_meta_highlight <- meta_var_list[meta_var_list %in% meta_data_highlight]
  
  # Abort if no meta_data_highlight found
  if (length(x = found_meta_highlight) == 0) {
    cli_abort(message = c("No 'meta_data_highlight' value(s) were not found.",
                          "i" = "The following {.code meta_data_highlight} variables were not found in {.field {good_meta_data_column}} and were omitted: {.field {bad_meta_highlight}}")
    )
  }
  
  # warn if some meta_data_highlight not found
  if (length(x = found_meta_highlight) != length(x = meta_data_highlight)) {
    cli_warn(message = c("Some {.code meta_data_highlight} value(s) were not found.",
                         "i" = "The following {.code meta_data_highlight} variables were not found in {.field {good_meta_data_column}} and were omitted: {.field {bad_meta_highlight}}")
    )
  }
  
  # Add raster check for scCustomize
  raster <- raster %||% (length(x = Cells(x = seurat_object)) > 2e5)
  
  # Change default ident and pull cells to highlight in plot
  Idents(object = seurat_object) <- good_meta_data_column
  
  cells_to_highlight <- CellsByIdentities(seurat_object, idents = found_meta_highlight)
  
  # set point size
  if (is.null(x = pt.size)) {
    pt.size <- AutoPointSize_scCustom(data = sum(lengths(x = cells_to_highlight)), raster = raster)
  }
  
  # Set colors
  # Adjust colors if needed when length(meta_data_highlight) > 1
  if (length(x = highlight_color) == 1 && length(x = found_meta_highlight) > 1) {
    highlight_color <- rep(x = highlight_color, length(x = found_meta_highlight))
    cli_inform(message = c("NOTE: Only one color provided to but {length(x = found_meta_highlight) `meta_data_highlight` variables were provided.}",
                           "i" = "Using the same color ({highlight_color[1]}) for all variables"))
  }
  
  # If NULL set using scCustomize_Palette
  if (is.null(x = highlight_color)) {
    highlight_color <- scCustomize_Palette(num_groups = length(x = cells_to_highlight), ggplot_default_colors = ggplot_default_colors)
  }
  
  # plot
  plot <- DimPlot_scCustom_man(seurat_object = seurat_object,
                           cells.highlight = cells_to_highlight,
                           cols.highlight = highlight_color,
                           colors_use = background_color,
                           sizes.highlight = pt.size,
                           pt.size = pt.size,
                           order = order,
                           raster = raster,
                           raster.dpi = raster.dpi,
                           split.by = split.by,
                           split_seurat = split_seurat,
                           label = label,
                           ...)
  
  # Update legend and return plot
  plot <- suppressMessages(plot & scale_color_manual(breaks = names(x = cells_to_highlight), values = c(highlight_color, background_color), na.value = background_color))
  
  # Aspect ratio changes
  if (!is.null(x = aspect_ratio)) {
    if (!is.numeric(x = aspect_ratio)) {
      cli_abort(message = "{.code aspect_ratio} must be a {.field numeric} value.")
    }
    plot <- plot & theme(aspect.ratio = aspect_ratio)
  }
  
  # Figure plot
  if (isTRUE(x = figure_plot)) {
    plot <- Figure_Plot(plot = plot)
  }
  
  return(plot)
}

# DimPlot_scCustom_man -------
# added order to param list
# changed DimPlot to DimPlot_man, added order in param
# removed raster = raster, raster.dpi = raster.dpi
DimPlot_scCustom_man <- function(
  seurat_object,
  colors_use = NULL,
  pt.size = NULL,
  reduction = NULL,
  group.by = NULL,
  split.by = NULL,
  split_seurat = FALSE,
  figure_plot = FALSE,
  aspect_ratio = NULL,
  shuffle = TRUE,
  seed = 1,
  label = NULL,
  label.size = 4,
  label.color = 'black',
  label.box = FALSE,
  dims = c(1, 2),
  repel = FALSE,
  raster = NULL,
  raster.dpi = c(512, 512),
  num_columns = NULL,
  ggplot_default_colors = FALSE,
  color_seed = 123,
  order = order,
  ...
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)
  
  # Change label if label.box
  if (isTRUE(x = label.box) && is.null(x = label)) {
    label <- TRUE
  }
  
  if (!is.null(x = split.by)) {
    split.by <- Meta_Present(object = seurat_object, meta_col_names = split.by, print_msg = FALSE, omit_warn = FALSE)[[1]]
  }
  
  # Add check for group.by before getting to colors
  if (length(x = group.by) > 1) {
    Meta_Present(object = seurat_object, meta_col_names = group.by, print_msg = FALSE)
  } else {
    if (!is.null(x = group.by) && group.by != "ident") {
      Meta_Present(object = seurat_object, meta_col_names = group.by, print_msg = FALSE)
    }
  }
  
  # Add one time split_seurat warning
  if (!is.null(x = split.by) && isFALSE(x = split_seurat) && getOption(x = 'scCustomize_warn_DimPlot_split_type', default = TRUE)) {
    cli_inform(c("",
                 "NOTE: {.field DimPlot_scCustom} returns split plots as layout of all plots each",
                 "with their own axes as opposed to Seurat which returns with shared x or y axis.",
                 "To return to Seurat behvaior set {.code split_seurat = TRUE}.",
                 "",
                 "-----This message will be shown once per session.-----"))
    options(scCustomize_warn_DimPlot_split_type = FALSE)
  }
  
  # Add raster check for scCustomize
  raster <- raster %||% (length(x = Cells(x = seurat_object)) > 2e5)
  
  label <- label %||% (is.null(x = group.by))
  
  # if split.by is null set split_seurat to TRUE
  if (is.null(x = split.by)) {
    split_seurat <- TRUE
  }
  
  # # figure_plot check
  # if (isTRUE(x = figure_plot) && isTRUE(x = split_seurat)) {
  #   cli_abort(message = "{.code figure_plot} can only be TRUE if {.code split_seurat = FALSE}.")
  # }
  
  # Set default color palette based on number of levels being plotted
  if (length(x = group.by) > 1) {
    all_length <- lapply(group.by, function(x) {
      num_var <- length(x = unique(x = seurat_object@meta.data[[x]]))
    })
    group_by_length <- max(unlist(x = all_length))
  } else {
    if (is.null(x = group.by)) {
      group_by_length <- length(x = unique(x = seurat_object@active.ident))
    } else {
      group_by_length <- length(x = unique(x = seurat_object@meta.data[[group.by]]))
    }
  }
  
  # Check colors use vs. ggplot2 color scale
  if (!is.null(x = colors_use) && isTRUE(x = ggplot_default_colors)) {
    cli_abort(message = "Cannot provide both custom palette to {.code colors_use} and specify {.code ggplot_default_colors = TRUE}.")
  }
  
  # set default plot colors
  if (is.null(x = colors_use)) {
    colors_use <- scCustomize_Palette(num_groups = group_by_length, ggplot_default_colors = ggplot_default_colors, color_seed = color_seed)
  }
  
  # Set uniform point size is pt.size = NULL (based on plot with most cells)
  if (is.null(x = pt.size) && !is.null(split.by)) {
    # cells per meta data
    cells_by_split <- data.frame(table(seurat_object@meta.data[, split.by]))
    # Identity with greatest number of cells
    max_cells <- max(cells_by_split$Freq)
    # modified version of the autopointsize function from Seurat
    pt.size <- AutoPointSize_scCustom(data = max_cells, raster = raster)
  }
  
  # set size otherwise
  pt.size <- pt.size %||% AutoPointSize_scCustom(data = seurat_object)
  
  # Plot
  if (is.null(x = split.by)) {
    plot <- DimPlot_man(object = seurat_object, cols = colors_use, pt.size = pt.size, reduction = reduction,
                        group.by = group.by, split.by = split.by, shuffle = shuffle, seed = seed, label = label, 
                        label.size = label.size, label.color = label.color, repel = repel, 
                        ncol = num_columns, dims = dims, label.box = label.box,
                        order = order, ...)
    if (isTRUE(x = figure_plot)) {
      
      # pull axis labels
      x_lab_reduc <- plot$labels$x
      y_lab_reduc <- plot$labels$y
      
      plot <- plot & NoAxes()
      
      axis_plot <- ggplot(data.frame(x= 100, y = 100), aes(x = .data[["x"]], y = .data[["y"]])) +
        geom_point() +
        xlim(c(0, 10)) + ylim(c(0, 10)) +
        theme_classic() +
        ylab(y_lab_reduc) + xlab(x_lab_reduc) +
        theme(plot.background = element_rect(fill = "transparent", colour = NA),
              panel.background = element_rect(fill = "transparent"),
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_line(
                arrow = arrow(angle = 15, length = unit(.5, "cm"), type = "closed")
              )
        )
      
      figure_layout <- c(
        area(t = 1, l = 2, b = 11, r = 11),
        area(t = 10, l = 1, b = 12, r = 2))
      
      plot_figure <- plot + axis_plot +
        plot_layout(design = figure_layout)
      
      # Aspect ratio changes
      if (!is.null(x = aspect_ratio)) {
        if (!is.numeric(x = aspect_ratio)) {
          cli_abort(message = "{.code aspect_ratio} must be a {.field numeric} value.")
        }
        plot_figure <- plot_figure & theme(aspect.ratio = aspect_ratio)
      }
      
      return(plot_figure)
    } else {
      # Aspect ratio changes
      if (!is.null(x = aspect_ratio)) {
        if (!is.numeric(x = aspect_ratio)) {
          cli_abort(message = "{.code aspect_ratio} must be a {.field numeric} value.")
        }
        plot <- plot & theme(aspect.ratio = aspect_ratio)
      }
      
      return(plot)
    }
    
  } else {
    if (isTRUE(x = split_seurat)) {
      # Plot Seurat Splitting
      plot <- DimPlot(object = seurat_object, cols = colors_use, pt.size = pt.size, reduction = reduction, group.by = group.by, split.by = split.by, shuffle = shuffle, seed = seed, label = label, label.size = label.size, label.color = label.color, repel = repel, raster = raster, raster.dpi = raster.dpi, ncol = num_columns, dims = dims, label.box = label.box, ...)
      if (isTRUE(x = figure_plot)) {
        
        # pull axis labels
        x_lab_reduc <- plot$labels$x
        y_lab_reduc <- plot$labels$y
        
        plot <- plot & NoAxes()
        
        axis_plot <- ggplot(data.frame(x= 100, y = 100), aes(x = .data[["x"]], y = .data[["y"]])) +
          geom_point() +
          xlim(c(0, 10)) + ylim(c(0, 10)) +
          theme_classic() +
          ylab(y_lab_reduc) + xlab(x_lab_reduc) +
          theme(plot.background = element_rect(fill = "transparent", colour = NA),
                panel.background = element_rect(fill = "transparent"),
                axis.text.x = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks = element_blank(),
                axis.line = element_line(
                  arrow = arrow(angle = 15, length = unit(.5, "cm"), type = "closed")
                )
          )
        
        figure_layout <- c(
          area(t = 1, l = 2, b = 11, r = 11),
          area(t = 10, l = 1, b = 12, r = 2))
        
        plot_figure <- plot + axis_plot +
          plot_layout(design = figure_layout)
        
        # Aspect ratio changes
        if (!is.null(x = aspect_ratio)) {
          if (!is.numeric(x = aspect_ratio)) {
            cli_abort(message = "{.code aspect_ratio} must be a {.field numeric} value.")
          }
          plot_figure <- plot_figure & theme(aspect.ratio = aspect_ratio)
        }
        
        return(plot_figure)
      } else {
        # Aspect ratio changes
        if (!is.null(x = aspect_ratio)) {
          if (!is.numeric(x = aspect_ratio)) {
            cli_abort(message = "{.code aspect_ratio} must be a {.field numeric} value.")
          }
          plot <- plot & theme(aspect.ratio = aspect_ratio)
        }
        
        return(plot)
      }
    } else {
      if (is.null(x = group.by)) {
        group_by_vars <- as.character(unique(x = seurat_object@active.ident))
      } else {
        group_by_vars <- as.character(unique(x = seurat_object@meta.data[[group.by]]))
      }
      # Extract reduction coordinates
      reduction <- reduction %||% DefaultDimReduc(object = seurat_object)
      all_cells <- Cells(x = seurat_object)
      reduc_coordinates <- Embeddings(object = seurat_object[[reduction]])[all_cells, dims]
      reduc_coordinates <- as.data.frame(x = reduc_coordinates)
      x_axis <- c(min(reduc_coordinates[, 1]),
                  max(reduc_coordinates[, 1]))
      y_axis <- c(min(reduc_coordinates[, 2]),
                  max(reduc_coordinates[, 2]))
      
      # Extract cell names per meta data list of values
      # Extract split.by list of values
      if (inherits(x = seurat_object@meta.data[, split.by], what = "factor")) {
        split_by_list <- as.character(x = levels(x = seurat_object@meta.data[, split.by]))
      } else {
        split_by_list <- as.character(x = unique(x = seurat_object@meta.data[, split.by]))
      }
      
      cell_names <- lapply(split_by_list, function(x) {
        row.names(x = seurat_object@meta.data)[which(seurat_object@meta.data[, split.by] == x)]})
      
      # Unify colors across plots
      if (is.null(x = group.by)) {
        levels_overall <- levels(x = Idents(object = seurat_object))
      } else {
        levels_overall <- levels(x = seurat_object@meta.data[[group.by]])
      }
      
      colors_overall <- colors_use
      
      names(x = colors_overall) <- levels_overall
      
      # plot
      plots <- lapply(1:length(x = split_by_list), function(x) {
        plot <- DimPlot(object = seurat_object, cells = cell_names[[x]], group.by = group.by, cols = colors_use, reduction = reduction, pt.size = pt.size, raster = raster, raster.dpi = raster.dpi, shuffle = shuffle, seed = seed, label = label, label.size = label.size, label.color = label.color, repel = repel, dims = dims, label.box = label.box, ...) +
          ggtitle(paste(split_by_list[[x]])) +
          theme(plot.title = element_text(hjust = 0.5),
                legend.position = "right") +
          xlim(x_axis) +
          ylim(y_axis)
        
        # Normalize the colors across all plots
        plot <- suppressMessages(plot + scale_color_manual(values = colors_overall, drop = FALSE))
        
        if (!is.null(x = group.by)) {
          plot <- plot + labs(color=group.by)
        } else {
          plot <- plot
        }
      })
      
      # Wrap Plots into single output
      plots <- wrap_plots(plots, ncol = num_columns) + plot_layout(guides = 'collect')
      
      # Aspect ratio changes
      if (!is.null(x = aspect_ratio)) {
        if (!is.numeric(x = aspect_ratio)) {
          cli_abort(message = "{.code aspect_ratio} must be a {.field numeric} value.")
        }
        plots <- plots & theme(aspect.ratio = aspect_ratio)
      }
      
      return(plots)
    }
  }
}



# DimPLot-----
#changed to singledimplot_man
DimPlot_man <- function(
  object,
  dims = c(1, 2),
  cells = NULL,
  cols = NULL,
  pt.size = NULL,
  reduction = NULL,
  group.by = NULL,
  split.by = NULL,
  shape.by = NULL,
  order = NULL,
  shuffle = FALSE,
  seed = 1,
  label = FALSE,
  label.size = 4,
  label.color = 'black',
  label.box = FALSE,
  repel = FALSE,
  cells.highlight = NULL,
  cols.highlight = '#DE2D26',
  sizes.highlight = 1,
  na.value = 'grey50',
  ncol = NULL,
  combine = TRUE
) {
  if (length(x = dims) != 2) {
    stop("'dims' must be a two-length vector")
  }
  reduction <- reduction %||% DefaultDimReduc(object = object)
  cells <- cells %||% colnames(x = object)
  if (isTRUE(x = shuffle)) {
    set.seed(seed = seed)
    cells <- sample(x = cells)
  }
  data <- Embeddings(object = object[[reduction]])[cells, dims]
  data <- as.data.frame(x = data)
  dims <- paste0(Key(object = object[[reduction]]), dims)
  object[['ident']] <- Idents(object = object)
  orig.groups <- group.by
  group.by <- group.by %||% 'ident'
  data[, group.by] <- object[[group.by]][cells, , drop = FALSE]
  for (group in group.by) {
    if (!is.factor(x = data[, group])) {
      data[, group] <- factor(x = data[, group])
    }
  }
  if (!is.null(x = shape.by)) {
    data[, shape.by] <- object[[shape.by, drop = TRUE]]
  }
  if (!is.null(x = split.by)) {
    data[, split.by] <- object[[split.by, drop = TRUE]]
  }
  plots <- lapply(
    X = group.by,
    FUN = function(x) {
      plot <- SingleDimPlot_man(
        data = data[, c(dims, x, split.by, shape.by)],
        dims = dims,
        col.by = x,
        cols = cols,
        pt.size = pt.size,
        shape.by = shape.by,
        order = order,
        label = FALSE,
        cells.highlight = cells.highlight,
        cols.highlight = cols.highlight,
        sizes.highlight = sizes.highlight,
        na.value = na.value
      )
      if (label) {
        plot <- LabelClusters(
          plot = plot,
          id = x,
          repel = repel,
          size = label.size,
          split.by = split.by,
          box = label.box,
          color = label.color
        )
      }
      if (!is.null(x = split.by)) {
        plot <- plot + FacetTheme() +
          facet_wrap(
            facets = vars(!!sym(x = split.by)),
            ncol = if (length(x = group.by) > 1 || is.null(x = ncol)) {
              length(x = unique(x = data[, split.by]))
            } else {
              ncol
            }
          )
      }
      return(plot)
    }
  )
  if (!is.null(x = split.by)) {
    ncol <- 1
  }
  if (combine) {
    plots <- wrap_plots(plots, ncol = orig.groups %iff% ncol)
  }
  return(plots)
}

# ----
# commented the order overide when cell selection is performed
SingleDimPlot_man <- function(
  data,
  dims,
  col.by = NULL,
  cols = NULL,
  pt.size = NULL,
  shape.by = NULL,
  alpha.by = NULL,
  order = NULL,
  label = FALSE,
  repel = FALSE,
  label.size = 4,
  cells.highlight = NULL,
  cols.highlight = '#DE2D26',
  sizes.highlight = 1,
  na.value = 'grey50'
) {
  pt.size <- pt.size %||% AutoPointSize(data = data)
  if (length(x = dims) != 2) {
    stop("'dims' must be a two-length vector")
  }
  if (!is.data.frame(x = data)) {
    data <- as.data.frame(x = data)
  }
  if (is.character(x = dims) && !all(dims %in% colnames(x = data))) {
    stop("Cannot find dimensions to plot in data")
  } else if (is.numeric(x = dims)) {
    dims <- colnames(x = data)[dims]
  }
  if (!is.null(x = cells.highlight)) {
    highlight.info <- SetHighlight(
      cells.highlight = cells.highlight,
      cells.all = rownames(x = data),
      sizes.highlight = sizes.highlight %||% pt.size,
      cols.highlight = cols.highlight,
      col.base = cols[1] %||% '#C3C3C3',
      pt.size = pt.size
    )
    # order <- highlight.info$plot.order
    data$highlight <- highlight.info$highlight
    col.by <- 'highlight'
    pt.size <- highlight.info$size
    cols <- highlight.info$color
  }
  if (!is.null(x = order) && !is.null(x = col.by)) {
    if (typeof(x = order) == "logical") {
      if (order) {
        data <- data[order(data[, col.by]), ]
      }
    } else {
      order <- rev(x = c(
        order,
        setdiff(x = unique(x = data[, col.by]), y = order)
      ))
      data[, col.by] <- factor(x = data[, col.by], levels = order)
      new.order <- order(x = data[, col.by])
      data <- data[new.order, ]
      if (length(x = pt.size) == length(x = new.order)) {
        pt.size <- pt.size[new.order]
      }
    }
  }
  if (!is.null(x = col.by) && !col.by %in% colnames(x = data)) {
    warning("Cannot find ", col.by, " in plotting data, not coloring plot")
    col.by <- NULL
  } else {
    # col.index <- grep(pattern = col.by, x = colnames(x = data), fixed = TRUE)
    col.index <- match(x = col.by, table = colnames(x = data))
    if (grepl(pattern = '^\\d', x = col.by)) {
      # Do something for numbers
      col.by <- paste0('x', col.by)
    } else if (grepl(pattern = '-', x = col.by)) {
      # Do something for dashes
      col.by <- gsub(pattern = '-', replacement = '.', x = col.by)
    }
    colnames(x = data)[col.index] <- col.by
  }
  if (!is.null(x = shape.by) && !shape.by %in% colnames(x = data)) {
    warning("Cannot find ", shape.by, " in plotting data, not shaping plot")
  }
  if (!is.null(x = alpha.by) && !alpha.by %in% colnames(x = data)) {
    warning(
      "Cannot find alpha variable ",
      alpha.by,
      " in data, setting to NULL",
      call. = FALSE,
      immediate. = TRUE
    )
    alpha.by <- NULL
  }
  plot <- ggplot(data = data) +
    geom_point(
      mapping = aes_string(
        x = dims[1],
        y = dims[2],
        color = paste0("`", col.by, "`"),
        shape = shape.by,
        alpha = alpha.by
      ),
      size = pt.size
    ) +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    labs(color = NULL)
  if (label && !is.null(x = col.by)) {
    plot <- LabelClusters(
      plot = plot,
      id = col.by,
      repel = repel,
      size = label.size
    )
  }
  if (!is.null(x = cols)) {
    if (length(x = cols) == 1 && (is.numeric(x = cols) || cols %in% rownames(x = brewer.pal.info))) {
      scale <- scale_color_brewer(palette = cols, na.value = na.value)
    } else if (length(x = cols) == 1 && (cols %in% c('alphabet', 'alphabet2', 'glasbey', 'polychrome', 'stepped'))) {
      colors <- DiscretePalette(length(unique(data[[col.by]])), palette = cols)
      scale <- scale_color_manual(values = colors, na.value = na.value)
    } else {
      scale <- scale_color_manual(values = cols, na.value = na.value)
    }
    plot <- plot + scale
  }
  plot <- plot + theme_cowplot()
  return(plot)
}

#SetHighlight ----
SetHighlight <- function(
  cells.highlight,
  cells.all,
  sizes.highlight,
  cols.highlight,
  col.base = 'black',
  pt.size = 1
) {
  if (is.character(x = cells.highlight)) {
    cells.highlight <- list(cells.highlight)
  } else if (is.data.frame(x = cells.highlight) || !is.list(x = cells.highlight)) {
    cells.highlight <- as.list(x = cells.highlight)
  }
  cells.highlight <- lapply(
    X = cells.highlight,
    FUN = function(cells) {
      cells.return <- if (is.character(x = cells)) {
        cells[cells %in% cells.all]
      } else {
        cells <- as.numeric(x = cells)
        cells <- cells[cells <= length(x = cells.all)]
        cells.all[cells]
      }
      return(cells.return)
    }
  )
  cells.highlight <- Filter(f = length, x = cells.highlight)
  names.highlight <- if (is.null(x = names(x = cells.highlight))) {
    paste0('Group_', 1L:length(x = cells.highlight))
  } else {
    names(x = cells.highlight)
  }
  sizes.highlight <- rep_len(
    x = sizes.highlight,
    length.out = length(x = cells.highlight)
  )
  cols.highlight <- c(
    col.base,
    rep_len(x = cols.highlight, length.out = length(x = cells.highlight))
  )
  size <- rep_len(x = pt.size, length.out = length(x = cells.all))
  highlight <- rep_len(x = NA_character_, length.out = length(x = cells.all))
  if (length(x = cells.highlight) > 0) {
    for (i in 1:length(x = cells.highlight)) {
      cells.check <- cells.highlight[[i]]
      index.check <- match(x = cells.check, cells.all)
      highlight[index.check] <- names.highlight[i]
      size[index.check] <- sizes.highlight[i]
    }
  }
  plot.order <- sort(x = unique(x = highlight), na.last = TRUE)
  plot.order[is.na(x = plot.order)] <- 'Unselected'
  highlight[is.na(x = highlight)] <- 'Unselected'
  highlight <- as.factor(x = highlight)
  return(list(
    plot.order = plot.order,
    highlight = highlight,
    size = size,
    color = cols.highlight
  ))
}
# is seurat ---
Is_Seurat <- function(
  seurat_object
) {
  if (!inherits(what = "Seurat", x = seurat_object)) {
    cli_abort(message = "{.code seurat_object} provided is not an object of class: Seurat.")
  }
}

# Meta Present ------
Meta_Present <- function(
  object,
  seurat_object = lifecycle::deprecated(),
  meta_col_names,
  print_msg = TRUE,
  omit_warn = TRUE,
  return_none = FALSE
) {
  # Check is slot is supplied
  if (lifecycle::is_present(seurat_object)) {
    lifecycle::deprecate_warn(when = "2.1.0",
                              what = "Meta_Present(seurat_object)",
                              with = "Meta_Present(object)",
                              details = c("!" = "Please adjust code now to prepare for full deprecation in v2.2.0.")
    )
    
  }
  
  # Set possible variables based on object type
  if (inherits(x = object, what = "Seurat")) {
    possible_features <- colnames(x = object@meta.data)
  }
  
  if (inherits(x = object, what = "liger")) {
    possible_features <- colnames(x = object@cell.data)
  }
  
  # If any features not found
  if (any(!meta_col_names %in% possible_features)) {
    bad_meta <- meta_col_names[!meta_col_names %in% possible_features]
    found_meta <- meta_col_names[meta_col_names %in% possible_features]
    
    if (isFALSE(return_none)) {
      if (length(x = found_meta) < 1) {
        cli_abort(message = c("No meta data columns found.",
                              "i" = "The following @meta.data columns were not found: {.field {glue_collapse_scCustom(input_string = bad_meta, and = TRUE)}}")
        )
      }
    }
    
    # Return message of features not found
    if (length(x = bad_meta) > 0 && isTRUE(x = omit_warn)) {
      cli_warn(message = c("The following @meta.data columns were omitted as they were not found:",
                           "i" = "{.field {glue_collapse_scCustom(input_string = bad_meta, and = TRUE)}}")
      )
    }
    
    # Return the found features omitting the not found ones.
    meta_list <- list(
      found_meta = found_meta,
      bad_meta = bad_meta
    )
    
    return(meta_list)
  }
  
  # Print all found message if TRUE
  if (isTRUE(x = print_msg)) {
    cli_inform(message = "All @meta.data columns present.")
  }
  
  # Return full input gene list.
  meta_list <- list(
    found_meta = meta_col_names,
    bad_meta = NULL
  )
  return(meta_list)
}

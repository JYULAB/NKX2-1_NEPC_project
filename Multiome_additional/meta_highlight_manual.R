SetHighlight_man <- function(
  cells.highlight,
  cells.all,
  sizes.highlight,
  cols.highlight,
  col.base = 'black',
  pt.size = 1
) {
  print("this ran")
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
  
  # highlight[is.na(x = highlight)] <- 'Unselected'
  # highlight.length <- length(unique(highlight))-1
  # highlight.order <- c(1:highlight.length, 'Unselected')
  # print(highlight.order)
  # # highlight <- factor(x = highlight, levels = highlight.order)
  # plot.order <- as.character(sort(x = unique(x = highlight), na.last = TRUE))
  # # plot.order <- order
  # print(plot.order)
  # plot.order[is.na(x = plot.order)] <- 'Unselected'
  # 
  return(list(
    plot.order = plot.order,
    highlight = highlight,
    size = size,
    color = cols.highlight
  ))
}

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
  print("in dimplot_sc")
  print(paste0("order in dimplt sc ", order))
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
    plot <- DimPlot_man(object = seurat_object, cols = colors_use, pt.size = pt.size, order = order,
                        reduction = reduction, group.by = group.by, split.by = split.by, 
                        shuffle = shuffle, seed = seed, label = label, label.size = label.size, 
                        label.color = label.color, repel = repel, ncol = num_columns, dims = dims, label.box = label.box, ...)
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
        plot <- DimPlot_man(object = seurat_object, cells = cell_names[[x]], group.by = group.by, cols = colors_use, reduction = reduction, pt.size = pt.size, raster = raster, raster.dpi = raster.dpi, shuffle = shuffle, seed = seed, label = label, label.size = label.size, label.color = label.color, repel = repel, dims = dims, label.box = label.box, ...) +
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

# highligh ----

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
  print(paste0("meta order", order))
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
  
  # return(list(cells_to_highlight, highlight_color))
  # plot
  print(paste0("order in meta ", order))
  
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



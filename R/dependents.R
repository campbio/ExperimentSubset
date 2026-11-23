# What happens to the subsets made from what goes away: a main assay, an
# assay of a subset or of an altExp, an altExp, or a subset (#122). Every
# path that can remove one takes dependents: "error" stops, "keep" makes
# each such subset store its own copy of what it saw, and "remove" removes
# them and those made from them.

.dependentsValues <- c("error", "keep", "remove")

# Stops unless dependents is one of .dependentsValues
.checkDependents <- function(dependents) {
  if (!.isOneString(dependents) || !dependents %in% .dependentsValues) {
    quoted <- sprintf("\"%s\"", .dependentsValues)
    stop("'dependents' must be ", quoted[1], ", ", quoted[2], " or ",
         quoted[3], ".")
  }
  invisible(dependents)
}

# The dependents argument in arglist, the arguments in ... of a setter,
# checked; "error" if not given
.dependentsArg <- function(arglist) {
  if (!"dependents" %in% names(arglist))
    return("error")
  .checkDependents(arglist[["dependents"]])
}

# assays<- on the main object, and so assay(x, i) <- NULL, which calls it:
# the subsets made from a removed assay follow dependents. A rename, such
# as from assayNames<-, which also calls it, isn't a removal (#126).
`.assays<-` <- function(x, withDimnames = TRUE, ..., value) {
  dependents <- .dependentsArg(list(...))
  out <- callNextMethod(x, withDimnames = withDimnames, ..., value = value)
  if (all(assayNames(x) %in% assayNames(out)) || .renamedOnly(x, out))
    return(out)
  .resolveDependents(x, out, dependents)
}

# TRUE if new has the assays of old, the same objects in the same order,
# under other names
.renamedOnly <- function(old, new) {
  a <- as.list(assays(old, withDimnames = FALSE))
  b <- as.list(assays(new, withDimnames = FALSE))
  length(a) == length(b) && all(mapply(identical, a, b))
}

# after, a change of before, with what dependents says done to the subsets
# that lost what they were made from (.lostParents()): an error naming the
# topmost of them, each keeping its own copy, or the topmost removed with
# those below them
.resolveDependents <- function(before, after, dependents = "error") {
  .checkDependents(dependents)
  if (subsetCount(after) == 0)
    return(after)
  lost <- .lostParents(before, after)
  if (length(lost) == 0)
    return(after)
  top <- lost[.topSubsets(before, names(lost))]
  if (dependents == "error")
    .stopLost(top)
  if (dependents == "remove")
    return(.removeLost(after, top, before))
  .keepLost(before, after, lost)
}

# The subsets of both before and after that, in after, no longer see what
# they were made from: the altExp they were made from with
# createSubset(altExp =), or their primary assay, where it was stored in
# before. A subset that now stores the assay itself, or sees it from a
# subset below where it was stored, hasn't lost it, and nor has one that
# stored it itself: it wasn't made from its own copy. A list named by subset
# of list(kind, assay, what): kind is "altExp" or "assay", assay the name
# of the lost assay, and what names what was lost, for messages. Chains
# are read in before: in after, a bare parentAssay naming a removed main
# assay could resolve to a subset that stores an assay of that name.
.lostParents <- function(before, after) {
  lost <- list()
  for (name in intersect(subsetNames(after), subsetNames(before))) {
    chain <- .subsetChain(before, name)
    oldRoot <- .rowRoot(before, name)
    root <- .renamedRoot(after, oldRoot)
    if (!is.null(root)
        && is.null(.rootShape(after, root$container, root$altExp))) {
      lost[[name]] <- list(kind = "altExp", assay = NULL,
                           what = paste0("altExp '", root$altExp, "'"))
      next
    }
    primary <- .primaryAssayName(before, name)
    owners <- .assayOwners(before, name, oldRoot, chain)
    if (!primary %in% names(owners))
      next
    owner <- owners[[primary]]
    if (owner == name)
      next
    holder <- owner
    altExp <- NULL
    if (owner == "") {
      holder <- .rootAssayHolders(before, oldRoot)[[primary]]
      altExp <- root$altExp
    }
    if (primary %in% .heldAssays(after, holder, altExp)
        || .heldCloser(after, name, root, chain, primary, owner))
      next
    lost[[name]] <- list(kind = "assay", assay = primary,
                         what = .assayWords(primary, holder, owner,
                                            oldRoot$altExp))
  }
  lost
}

# Names of the assays that holder ("" for the main object, or a subset)
# stores in x: its own, or with altExp those of its altExp of that name
.heldAssays <- function(x, holder, altExp = NULL) {
  if (holder == "") {
    object <- x
  } else {
    if (!holder %in% subsetNames(x))
      return(character(0))
    object <- .internalAssay(.subsets(x)[[holder]])
  }
  if (is.null(altExp))
    return(assayNames(object))
  if (!altExp %in% .storedAltExpNames(object))
    return(character(0))
  assayNames(.storedAltExp(object, altExp))
}

# The row root of a subset in before (from .rowRoot()), with the name its
# altExp has in after, which may have renamed it; NULL for NULL
.renamedRoot <- function(after, root) {
  if (is.null(root))
    return(NULL)
  root$altExp <- .rootAltExp(.subsets(after)[[root$start]])
  root
}

# TRUE if subset name, with row root root and chain (from before) in x,
# sees assayName from a subset of its chain below owner, where it saw it
# before ("" for the root): a copy stored since, not a loss
.heldCloser <- function(x, name, root, chain, assayName, owner) {
  owners <- .assayOwners(x, name, root, chain)
  if (!assayName %in% names(owners) || owners[[assayName]] == "")
    return(FALSE)
  if (owner == "")
    return(TRUE)
  chain <- .rootChain(chain, root)
  isTRUE(match(owners[[assayName]], chain) > match(owner, chain))
}

# How a message names a lost assay: stored in holder, where owner is ""
# when it was linked from the root of the rows, and altExp the name of
# the altExp the rows come from, if any
.assayWords <- function(assayName, holder, owner, altExp) {
  if (owner != "")
    return(paste0("assay '", assayName, "' of subset '", holder, "'"))
  if (!is.null(altExp))
    assayName <- paste0(altExp, "/", assayName)
  words <- paste0("assay '", assayName, "'")
  if (holder != "")
    words <- paste0(words, " of subset '", holder, "'")
  words
}

# The names in names that no subset up their chain in x is also in
.topSubsets <- function(x, names) {
  Filter(function(name) {
    chain <- .subsetChain(x, name)
    !any(chain[-length(chain)] %in% names)
  }, names)
}

# "subset 'a'" or "subsets 'a' and 'b'"
.subsetWords <- function(names) {
  paste0(if (length(names) == 1) "subset " else "subsets ",
         .andList(sprintf("'%s'", names)))
}

# The end of the error for dependents = "error", for n subsets; without
# keep, only "remove" is offered
.dependentsHint <- function(n, keep = TRUE) {
  it <- if (n == 1) c("it", "its") else c("them", "their")
  remove <- paste0("dependents = \"remove\" to remove ", it[1],
                   " and the subsets made from ", it[1], ".")
  if (!keep)
    return(paste0(" Use ", remove))
  paste0(" Use dependents = \"keep\" to keep ", it[1], ", with ", it[2],
         " own copy of that data, or ", remove)
}

# Stops, naming each subset in lost (from .lostParents()) and what it was
# made from
.stopLost <- function(lost) {
  whats <- vapply(lost, `[[`, character(1), "what")
  sentences <- vapply(unique(whats), function(what) {
    subsets <- names(lost)[whats == what]
    kind <- lost[[subsets[1]]]$kind
    if (kind != "altExp")
      what <- paste0(toupper(substr(what, 1, 1)), substring(what, 2))
    paste0(what, " is the parent", if (kind == "assay") " assay", " of ",
           .subsetWords(subsets), ".")
  }, character(1))
  altExps <- vapply(lost, function(l) l$kind == "altExp", logical(1))
  stop(paste(sentences, collapse = " "),
       .dependentsHint(length(lost), keep = !any(altExps)))
}

# x without the subsets in lost (from .lostParents()) and those below them
# in chains (the object their chains are read in), with a warning for each
# thing lost that names them
.removeLost <- function(x, lost, chains = x) {
  whats <- vapply(lost, `[[`, character(1), "what")
  for (what in unique(whats)) {
    top <- names(lost)[whats == what]
    gone <- c(top, unlist(lapply(top, function(name) {
      .subsetsBelow(chains, name)
    })))
    gone <- subsetNames(x)[subsetNames(x) %in% gone]
    x <- .withoutSubsets(x, gone)
    warning("Removed ", .subsetWords(gone), ": the ", what,
            if (length(gone) == 1) " it was" else " they were",
            " made from was removed.", call. = FALSE)
  }
  x
}

# x without the subsets in names
.withoutSubsets <- function(x, names) {
  .subsets(x) <- .subsets(x)[!subsetNames(x) %in% names]
  x
}

# after, with each subset in lost (from .lostParents()) storing as its own
# what it saw of its lost assay in before, with a message; parents first,
# by their chains in before, so one that sees its parent's copy needs
# none. Stops, before anything
# is stored, if an altExp they were made from was removed (#127).
.keepLost <- function(before, after, lost) {
  if (any(vapply(lost, function(l) l$kind == "altExp", logical(1))))
    stop("dependents = \"keep\" isn't supported yet when an altExp is ",
         "removed (#127). Use dependents = \"remove\".")
  chains <- lapply(names(lost), function(name) .subsetChain(before, name))
  names(chains) <- names(lost)
  kept <- character(0)
  for (name in names(lost)[order(lengths(chains))]) {
    assayName <- lost[[name]]$assay
    chain <- chains[[name]]
    root <- .renamedRoot(after, .rowRoot(before, name))
    if (any(chain[-length(chain)] %in% kept)
        && assayName %in% names(.assayOwners(after, name, root, chain)))
      next
    after <- .storeOwnAssay(after, name, assayName,
                            .subsetAssayView(before, name, assayName,
                                             withDimnames = FALSE))
    message("Subset '", name, "' now stores its own assay '", assayName,
            "'.")
    kept <- c(kept, name)
  }
  after
}

# x with value stored as subset name's own assay assayName
.storeOwnAssay <- function(x, name, assayName, value) {
  internal <- .internalAssay(.subsets(x)[[name]])
  assay(internal, assayName, withDimnames = FALSE) <- value
  .internalAssay(.subsets(x)[[name]]) <- internal
  x
}

# removeSubset(): x without the subsets in subsetName, with dependents
# saying what happens to the subsets made from them
.removeSubset <- function(x, subsetName, dependents = "error") {
  .checkDependents(dependents)
  if (!is.character(subsetName) || length(subsetName) == 0
      || anyNA(subsetName))
    stop("'subsetName' must be one or more subset names.")
  for (name in subsetName)
    .isSubset(x, name)
  x <- .qualifyParents(x)
  named <- unique(subsetName)
  below <- unique(unlist(lapply(named, function(name) {
    .subsetsBelow(x, name)
  })))
  # the subsets made directly from a removed one that stay
  children <- setdiff(below, named)
  children <- children[vapply(children, function(name) {
    .parentSubset(x, name) %in% named
  }, logical(1))]
  if (length(children) == 0)
    return(.withoutSubsets(x, named))
  if (dependents == "keep")
    return(.removeKeeping(x, named))
  lost <- lapply(children, function(name) {
    list(kind = "subset", assay = NULL,
         what = paste0("subset '", .parentSubset(x, name), "'"))
  })
  names(lost) <- children
  if (dependents == "error")
    .stopLost(lost)
  .withoutSubsets(.removeLost(x, lost), named)
}

# The subset subsetName was made from, or NULL for the main object
.parentSubset <- function(x, subsetName) {
  chain <- .subsetChain(x, subsetName)
  if (length(chain) == 1)
    return(NULL)
  chain[length(chain) - 1]
}

# x without the subsets in named, each removed in turn from the top down,
# after the subsets made directly from it are made from its parent instead
# (.reRoot()) and store what they saw from it (.keepSeen()). Subsets that
# are removed later say nothing about what they stored.
.removeKeeping <- function(x, named) {
  depth <- vapply(named, function(name) length(.subsetChain(x, name)),
                  integer(1))
  for (removed in named[order(depth)]) {
    children <- Filter(function(name) {
      identical(.parentSubset(x, name), removed)
    }, subsetNames(x))
    before <- x
    for (child in children)
      x <- .reRoot(before, x, child, removed)
    x <- .withoutSubsets(x, removed)
    for (child in children)
      x <- .keepSeen(before, x, child, removed,
                     quiet = child %in% named)
  }
  x
}

# x with subset child, made directly from subset removed, made from
# removed's parent instead: its indices composed through removed's, and
# its parentAssay naming its primary assay one level up. A child made
# with createSubset(altExp =) from removed's view of an altExp keeps its
# rows if removed's parent sees that altExp with the same rows; one below
# a subset made with createSubset(altExp =) that is removed becomes one.
.reRoot <- function(before, x, child, removed) {
  subset <- .subsets(x)[[child]]
  parent <- .subsets(x)[[removed]]
  up <- .parentSubset(before, removed)
  cols <- .colIndices(parent)[.colIndices(subset)]
  rows <- .rowIndices(parent)[.rowIndices(subset)]
  altExp <- .rootAltExp(subset)
  if (!is.null(altExp)) {
    # its rows are positions in that altExp's features
    rows <- .rowIndices(subset)
    .checkSameRoot(before, child, removed, up, altExp)
  } else {
    altExp <- .rootAltExp(parent)
  }
  primary <- .primaryAssayName(before, child)
  parentAssay <- if (is.null(up)) primary else paste0(up, "/", primary)
  # made from what removed was made from, linked: made from that too
  if (is.null(.rootAltExp(subset))
      && .madeFromLinked(before, removed, primary))
    parentAssay <- .parentAssay(parent)
  .subsets(x)[[child]] <- .copySubset(subset, rowIndices = rows,
                                      colIndices = cols,
                                      parentAssay = parentAssay,
                                      altExpName = altExp)
  x
}

# TRUE if primary is the primary assay of subset removed, which it sees
# linked rather than stores itself
.madeFromLinked <- function(x, removed, primary) {
  if (is.null(.parentAssay(.subsets(x)[[removed]]))
      || !identical(.primaryAssayName(x, removed), primary))
    return(FALSE)
  owners <- .assayOwners(x, removed)
  if (!primary %in% names(owners))
    return(FALSE)
  owner <- owners[[primary]]
  if (owner == "")
    owner <- .rootAssayHolders(x, .rowRoot(x, removed))[[primary]]
  owner != removed
}

# Stops unless up (a subset, or NULL for the main object) sees the altExp
# called altExp with the rows that removed sees, for child, made from
# removed's view of it
.checkSameRoot <- function(x, child, removed, up, altExp) {
  old <- .rootShape(x, removed, altExp)
  new <- .rootShape(x, up, altExp)
  if (is.null(new) || !.sameAltExpShape(old, new))
    stop("Subset '", child, "' can't be kept: altExp '", altExp,
         "' of subset '", removed, "' has other rows than that of its ",
         "parent. Use dependents = \"remove\".")
  invisible(NULL)
}

# x with subset child, re-rooted from subset removed (.reRoot()), storing
# as its own what it saw in before from removed's own storage: assays,
# altExps, spatial coordinates and trees, with a message unless quiet
.keepSeen <- function(before, x, child, removed, quiet = FALSE) {
  say <- if (quiet) function(...) invisible(NULL) else message
  # assays first, so a parentAssay that names one of them finds it
  stored <- character(0)
  owners <- .assayOwners(before, child)
  rootHolders <- .rootAssayHolders(before, .rowRoot(before, child))
  for (name in names(owners)) {
    owner <- owners[[name]]
    from <- if (owner == "") rootHolders[[name]] else owner
    if (from != removed)
      next
    x <- .storeOwnAssay(x, child, name,
                        .subsetAssayView(before, child, name,
                                         withDimnames = FALSE))
    stored <- c(stored, name)
  }
  altExps <- .keepSeenAltExps(before, x, child, removed)
  x <- altExps$x
  stored <- c(stored, altExps$stored)
  if (length(stored) > 0)
    say("Subset '", child, "' now stores its own ",
        if (length(stored) == 1) "assay " else "assays ",
        .andList(sprintf("'%s'", stored)), ".")
  if (identical(.spatialCoordsOwner(before, child), removed)) {
    coords <- .subsetSpatialCoords(before, child)
    rownames(coords) <- NULL
    .subsetSpatialCoords(x, child) <- coords
    say("Subset '", child, "' now stores its own spatialCoords.")
  }
  for (dim in c("row", "col")) {
    if (dim == "row" && !is.null(.rowRoot(before, child)))
      next
    if (identical(.treeOwner(before, child, dim), removed)) {
      .subsetTree(x, child, dim) <- .subsetTree(before, child, dim)
      say("Subset '", child, "' now stores its own ", .treeDim(dim)$name,
          ".")
    }
  }
  x
}

# x with subset child storing what it saw in before of the altExps that
# removed stored, so it sees them as before: a list of that x and stored,
# the names of the altExp assays it now stores, as "altExp/assay"
.keepSeenAltExps <- function(before, x, child, removed) {
  unchanged <- list(x = x, stored = character(0))
  if (!.holdsSce(before, child))
    return(unchanged)
  visible <- .altExpHolders(before, child)
  seen <- vapply(visible$holders, function(holders) removed %in% holders,
                 logical(1))
  if (!any(seen))
    return(unchanged)
  ownBefore <- .ownAltExpAssays(x, child)
  views <- SingleCellExperiment::altExps(.withSubsetAltExps(before, child),
                                         withDimnames = FALSE)
  x <- .storeAltExps(x, child, views, unlink = TRUE, quiet = TRUE)
  list(x = x, stored = setdiff(.ownAltExpAssays(x, child), ownBefore))
}

# The assays of the altExps that subset name stores itself, as
# "altExp/assay"
.ownAltExpAssays <- function(x, name) {
  internal <- .internalAssay(.subsets(x)[[name]])
  out <- character(0)
  for (altExp in .storedAltExpNames(internal))
    out <- c(out, paste0(altExp, "/",
                         assayNames(.storedAltExp(internal, altExp))))
  out
}

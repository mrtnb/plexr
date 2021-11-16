
imbalance_penalty_function1 <- function(x) exp(7*(x^4)) + (20*x) -1

score_plex <- function(plex, ...){
  norm_imbalance <- abs(channel_A_freq(plex, ...) - 0.5) / 0.5
  return(-1 * sum(imbalance_penalty_function1(norm_imbalance)))
}

channel_AB_counts <- function(plex, channel_A=c('A','C'), channel_B=c('G', 'T')){
  apply(plex, MARGIN=2, function(cycle){
    channel_bias <- c('A'=sum(cycle %in% channel_A), 'B'=sum(cycle %in% channel_B))
    return(channel_bias)
  })
}

channel_A_freq <- function(plex, ...){
  apply(channel_AB_counts(plex, ...), MARGIN=2, function(x) x['A']/sum(x))
}


score_plex_list <- function(plex_list_obj, ...){
  plex_list <- plex_list_obj$pl
  s<-sum(sapply(plex_list, score_plex))
  if(is.nan(s)){
    plex_list_summary(plex_list)
    stop()
  }
  return(s)
  #return(weighted.mean(x=sapply(plex_list, score_plex), w=sapply(plex_list, length)))
}

sample.safe <- function(x, ...){
  ret <- NULL
  if(length(x)==1){
    ret <- x
  }else{
    ret <- sample(x, ...)
  }
  return(ret)
}

mutate_plex_list <- function(plex_list_obj, prob_split=0.5, prob_merge=0.5, prob_move=0.5, min_plex_size=2, max_plex_size=12){
  plex_list <- plex_list_obj$pl
  # message('Mut')
  # message('prob_split: ', prob_split)
  # message('prob_merge: ', prob_merge)
  # message('prob_move: ', prob_move)
  # message('min_plex_size: ', min_plex_size)
  # message('max_plex_size: ', max_plex_size)
  rnames_start <- sort(rownames(do.call(rbind, plex_list)))
  n_list = length(plex_list)
  plex_sizes <- sapply(plex_list, nrow)
  
  if(runif(1)<prob_split & any(plex_sizes>min_plex_size)){
    #message('Split')
    # TODO plex_sizes > 2*min_plex_size ?
    split_idx <- sample.safe(which(plex_sizes>min_plex_size), 1)
    plex_split <- plex_list[[split_idx]]
    stopifnot(nrow(plex_split)>min_plex_size)
    # TODO split_plex() does not obey min split size.....
    plex_list <- c(plex_list[-split_idx], split_plex(plex_split))
    n_list = length(plex_list)
    plex_sizes <- sapply(plex_list, nrow)  
  }

  if(runif(1)<prob_merge & length(plex_list)>1){
    #message('Merge')
    # collect plex index combinations of potential merges
    index_comb<-unlist(lapply(1:(n_list-1), function(x){
      lapply((x+1):n_list, function(y){
        c(x,y)
      })
    }), recursive = F)
    # Determine potential merge sizes
    merge_sizes <- sapply(index_comb, function(x){
      nrow(plex_list[[x[1]]]) + nrow(plex_list[[x[2]]])
    })
    if(any(merge_sizes <= max_plex_size)){
      merge_idxs <- index_comb[[sample.safe(which(merge_sizes <= max_plex_size), 1)]]
      merged_plex <- rbind(plex_list[[merge_idxs[1]]], plex_list[[merge_idxs[2]]])
      plex_list <- c(plex_list[-merge_idxs], list(merged_plex))
      n_list = length(plex_list)
      plex_sizes <- sapply(plex_list, nrow)  
    }
  }
  
  if(runif(1)<prob_move && length(plex_list)>1 && any(plex_sizes>min_plex_size) && any(plex_sizes<max_plex_size) && 
     sum(plex_sizes>min_plex_size|plex_sizes<max_plex_size)>=2){
    #message('Move')
    # Collect possible index combinations for moves
    index_comb<-unlist(lapply(1:n_list, function(x){
      lapply(1:n_list, function(y){
        c(x,y)
      })
    }), recursive = F)    
    # Determine valid move combinations
    valid_move_comb <- sapply(index_comb, function(x){
      x[1]!=x[2] && plex_sizes[x[1]]>min_plex_size && plex_sizes[x[2]]<max_plex_size
    })
    move_idx_comb <- sample.safe(index_comb[valid_move_comb], 1)[[1]]
    from_plex_idx <- move_idx_comb[1]
    to_plex_idx <- move_idx_comb[2]
    #to_plex_idx <- sample.safe(setdiff(1:length(plex_list), from_plex_idx),1)
    from_plex_primer_idx <- sample.safe(1:nrow(plex_list[[from_plex_idx]]), 1)
    primer <- plex_list[[from_plex_idx]][from_plex_primer_idx,,drop=F]
    new_from_plex <- plex_list[[from_plex_idx]][-from_plex_primer_idx,,drop=F]
    new_to_plex <- rbind(plex_list[[to_plex_idx]], primer)
    plex_list[[from_plex_idx]] <- new_from_plex
    plex_list[[to_plex_idx]] <- new_to_plex
    plex_sizes <- sapply(plex_list, nrow)  
  }
  rnames_end <- sort(rownames(do.call(rbind, plex_list)))
  stopifnot(all(rnames_start==rnames_end))
  return(PlexList(plex_list))
}

split_plex <- function(plex_split){
  plex_split_n <- nrow(plex_split)
  stopifnot(plex_split_n>1)
  rnd_split <- rep(T,plex_split_n)
  while(sum(rnd_split)==length(rnd_split) || sum(!rnd_split)==length(rnd_split)){
    rnd_split <- sample(c(T,F), plex_split_n, replace=T)
  }
  plex_list <- c(list(plex_split[rnd_split,,drop=FALSE]), 
                 list(plex_split[!rnd_split,,drop=FALSE]))
  return(plex_list)
}

create_random_split_plexlist <- function(proto_plexlist_obj, min_plex_size, max_plex_size, try_limit=10000){
  proto_plexlist <- proto_plexlist_obj$pl
  single_plex <- do.call(rbind, proto_plexlist)
  n_plex <- nrow(single_plex)
  plex_splits <- c(0)
  while(sum(plex_splits) != n_plex){
    next_split <- sample.safe(min_plex_size:max_plex_size, 1)
    this_try <- 1
    while(this_try <= try_limit && sum(plex_splits) + next_split != n_plex && 
          (sum(plex_splits) + next_split > n_plex || sum(plex_splits) + next_split > n_plex - min_plex_size)){
      next_split <- sample.safe(min_plex_size:max_plex_size, 1)
      this_try<-this_try+1
    }
    if(this_try >= try_limit){
      stop('Failed creating plex list. Exeeded ', try_limit, 
           ' iterations. min_plex_size=', min_plex_size, ' max_plex_size=', 
           max_plex_size, ' n_plex=', n_plex, ' plex_splits=', paste(plex_splits, collapse=','))
    }
    plex_splits <- c(plex_splits, next_split)
  }
  plex_splits <- plex_splits[plex_splits>0]
  new_plex_list <- split.data.frame(single_plex, rep(1:length(plex_splits), plex_splits))
  names(new_plex_list) <- NULL
  return(PlexList(new_plex_list))
}

bootstrap_population <- function(plex_list_obj, pop_size, min_plex_size, max_plex_size){
  pop <- list()
  #pop[[1]] <- mutate_plex_list(list(do.call(rbind, plex_list)), min_plex_size=min_plex_size)
  pop[[1]] <- create_random_split_plexlist(plex_list_obj, min_plex_size, max_plex_size)
  for(i in 2:pop_size){
    prob_split <- ((2*pop_size)-i)/(2*pop_size)
    pop[[i]] <- mutate_plex_list(pop[[i-1]], prob_split = prob_split, prob_merge = 1-prob_split, 
                                 min_plex_size=min_plex_size, max_plex_size=max_plex_size)
  }
  return(pop)
}

evolve_plex_list <- function(plex_list_obj, n_iter=10, n_gen=100, pop_size=100, death_rate=0.1, min_plex_size=2, max_plex_size=12, ...){
  winners <- lapply(1:n_iter, function(iter){
    #plex_list_pop <- c(list(plex_list), lapply(1:(pop_size-1), function(x) mutate_plex_list(plex_list, ...)))
    plex_list_pop <- bootstrap_population(plex_list_obj, pop_size, min_plex_size=min_plex_size, max_plex_size=max_plex_size)
    for(gen in 1:n_gen){
      plex_list_scores <- sapply(plex_list_pop, score_plex_list)
      message(sprintf('Iter %s, gen %s, max.score = %s', iter, gen, max(plex_list_scores)))
      death_idxs <- order(plex_list_scores)[1:(round(death_rate*length(plex_list_scores)))]
      death_sel <- 1:length(plex_list_scores) %in% death_idxs
      n_death <- length(death_idxs)
      prolif_sel <- sample.safe(which(!death_sel), n_death, replace=F)
      prolif_survivors <- plex_list_pop[prolif_sel]
      new_pop <- lapply(prolif_survivors, mutate_plex_list, min_plex_size=min_plex_size, max_plex_size=max_plex_size, ...)
      plex_list_pop <- c(plex_list_pop[!death_sel], new_pop)
    }
    iter_plex_list_scores <- sapply(plex_list_pop, score_plex_list)
    iter_winner <- plex_list_pop[[which.max(iter_plex_list_scores)]]
    return(iter_winner)
  })
  plex_list_scores <- sapply(winners, score_plex_list)
  winner <- winners[[which.max(plex_list_scores)]]
  message('Best plex list has score: ', score_plex_list(winner))
  return(name_sort_plex_list(winner))
}

name_sort_plex_list <- function(plex_list_obj){
  plex_list_obj$pl <- lapply(plex_list_obj$pl, function(x) x[order(rownames(x)),])
  return(plex_list_obj)
}

reduce_to_plex <- function(plexlike){
  reduced <- NULL
  if(class(plexlike)=='PlexList'){
    reduced <- do.call(rbind, plexlike$pl)
  }else if(class(plexlike)=='list'){
    reduced <- do.call(rbind, plexlike)
  } else {
    reduced <- plexlike
  }
  return(reduced)
}

find_plex_friends <- function(target_plex, pool_plex, n_friends, max_tests=100000){
  target_plex <- reduce_to_plex(target_plex)
  pool_plex <- reduce_to_plex(pool_plex)
  n_tests <- choose(nrow(pool_plex), n_friends)
  if(n_tests <= max_tests){
    pool_plex_combs <- combn(nrow(pool_plex), n_friends, simplify=F, FUN=function(test_friends){
      pool_plex[test_friends,,drop=F]
    })
    extended_target_plex <- lapply(pool_plex_combs, function(pool_plex_select){
      rbind(target_plex, pool_plex_select)
    })
    scores <- sapply(extended_target_plex, score_plex)
    max.score <- max(scores)
    list(friends=pool_plex_combs[scores==max.score], plexes=extended_target_plex[scores==max.score], score=max.score)
  } else {
    stop(sprintf('Number required tests (N=%s) exeeds max_tests=%s', n_tests, max_tests))
  }
}

read_plex_list <- function(filename){
  file_str <- readChar(filename, file.info(filename)$size)
  plexes_str_list <- strsplit(file_str,"\n\n")[[1]]
  plex_list <- lapply(plexes_str_list, function(plex_str){
    plex_lines <- strsplit(plex_str,"\n")[[1]]
    plex_struct <- sapply(plex_lines, function(primer_str){
      cols <- strsplit(primer_str, '\\s+')[[1]]
      if(length(cols)>1){
        primer_str <- cols[2]
      }
      return(strsplit(primer_str, ""))
    })
    names(plex_struct) <- sapply(strsplit(plex_lines, '\\s+'), function(x) x[1])
    return(do.call(rbind, plex_struct))
  })
  return(PlexList(plex_list))
}

write_plex_list <- function(plex_list_obj, filename){
  plex_list <- plex_list_obj$pl
  file_con <- file(filename)
  first = TRUE
  line_buffer <- character()
  for (i in 1:length(plex_list)){
    plex <- plex_list[[i]]
    if(i>1){
      line_buffer <- c(line_buffer, "")
    }
    seqs <- apply(plex, 1, function(x) paste(x, collapse=''))
    if(!is.null(names(seqs))){
      line_buffer <- c(line_buffer, paste(names(seqs), seqs, sep="\t"))
    }else{
      line_buffer <- c(line_buffer, as.character(seqs))
    }
  }
  text<-paste(line_buffer, collapse="\n")
  writeLines(text, file_con, sep='')
  close(file_con)
}

PlexList <- function(pl){
  return(structure(list(pl=pl), class='PlexList'))
}

print.PlexList <- function(plex_list_obj){
  plex_list <- plex_list_obj$pl
  n_plex <- length(plex_list)
  cat(sprintf('List of %s plexes', n_plex), '\n')
  for(i in 1:n_plex){
    plex <- plex_list[[i]]
    cat(sprintf('%s. / %s-plex', i, nrow(plex)), '\n')
    print(plex)
    channel_counts_str <- paste(sapply(apply(channel_AB_counts(plex), 2, paste, collapse=','), 
                                       function(x) sprintf('(%s)', x)), collapse=' ')
    cat('-> Channel A/B counts :', channel_counts_str, '\n')
    cat('-> Channel A freq.    :', paste(round(channel_A_freq(plex), 2), collapse=', '), '\n')
    cat('-> Score         :', score_plex(plex), '\n')
    cat('\n')
  }
  cat('=> Overall score of plex list:', score_plex_list(plex_list_obj), '\n')
}

# test_plex1 <- matrix(c('A','C','G',
#                  'A','C','G',
#                  'A','C','G',
#                  'A','C','G'), nrow=4, byrow = T)
# 
# test_plex2 <- matrix(c('G','T','T',
#                   'G','T','A',
#                   'G','T','T',
#                   'G','T','T'), nrow=4, byrow = T)
# 
# test_plex_list <- list(test_plex1, test_plex2)
# plex_list_summary(test_plex_list)
# 
# test_winner <-evolve_plex_list(test_plex_list)
# score_plex_list(test_winner)
# 



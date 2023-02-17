

res <- addAttributes(data_tree = b_tree, df = aer)
res$Do(asrUpstream, traversal = 'post-order')
res$Do(inhDownstream, traversal = 'pre-order')


# jj <- res$p__1224$c__1236$o__91347$f__543$g__1903434
# # jj$aerobic__Score <- NULL
# # jj$facultatively_anaerobic__Score <- NULL

# jj <- res$p__57723$c__204432$o__204433$f__204434
# jj$Do(.my_fun, traversal = 'post-order')

print(
    # res$p__57723$c__204432$o__204433$f__204434,
    # res$p__1224$c__1236$o__91347$f__543,
    # res$p__1224$c__1236$o__91347$f__543$g__561,
    # jj,
    res,
    'aerobic__Score', 'facultatively_anaerobic__Score',
    'anaerobic__Score', 'microaerophilic__Score',
    'obligately_aerobic__Score', 'obligately_anaerobic__Score',
    limit = 1000
) |> View()

# res$p__1224$c__1236$o__91347$f__543$g__1903434
# res$p__57723$c__204432$o__204433$f__204434

print(
    # res$p__57723$c__204432$o__204433$f__204434,
    # res$p__1224$c__1236$o__91347$f__543,
    # res$p__1224$c__1236$o__91347$f__543$g__561,
    # jj,
    res,
    'aerobic__Score', 'aerobic__Evidence',
    'facultatively_anaerobic__Score', 'facultatively_anaerobic__Evidence',
    'anaerobic__Score', 'anaerobic__Evidence',
    'microaerophilic__Score', 'microaerophilic__Evidence',
    'obligately_aerobic__Score', 'obligately_aerobic__Evidence',
    'obligately_anaerobic__Score', 'obligately_anaerobic__Evidence',
    limit = 1000
) |> View()

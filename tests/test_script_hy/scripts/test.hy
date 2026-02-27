(require hyrule [-> ->>])

(defn is-odd? [n] (!= (% n 2) 0))

(setv result
      (->> (get snakemake.input 0)
           open
           .readlines
           (map int)
           (filter is-odd?)
           sum))

(print result
       :file
       (-> (get snakemake.output 0)
           (open "w")))

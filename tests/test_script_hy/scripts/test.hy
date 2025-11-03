(require hyrule [-> ->>])

(defn is-odd? [n] (!= (% n 2) 0))

(setv file (open "test.in"))
(with [file
       (-> (get snakemake.input 0)
           open)]
      (setv result
            (->> file
                 .readlines
                 (map int)
                 (filter is-odd?)
                 sum))
      (print result
             :file
             (-> (get snakemake.output 0)
                 (open "w"))))

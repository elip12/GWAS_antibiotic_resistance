CREATE TABLE IF NOT EXISTS kmer (
    seq         VARCHAR(30) NOT NULL,
    sample_id   VARCHAR(24) NOT NULL,
    contig      INT NOT NULL,
    index       INT NOT NULL,
    PRIMARY KEY(seq, sample_id, contig, index)
);

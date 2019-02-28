-- -- CREATE DOMAIN varchar_id AS
-- --     VARCHAR(24) NOT NULL

-- -- CREATE DOMAIN serial_id AS
-- --     SERIAL NOT NULL

-- -- CREATE TABLE IF NOT EXISTS sample (
-- --     sample_id           varchar_id UNIQUE     PRIMARY KEY
-- -- );

-- -- CREATE TABLE IF NOT EXISTS contig (
-- --     contig_id           serial_id UNIQUE          PRIMARY KEY,
-- --     contig_sequence     VARCHAR(10000) NOT NULL
-- -- );

-- -- CREATE TABLE IF NOT EXISTS sample_contig (
-- --     sample_contig_id    serial_id UNIQUE          PRIMARY KEY,
-- --     sample_id           VARCHAR(24) NOT NULL REFERENCES sample,
-- --     contig_id           serial_id REFERENCES contig
-- -- );

-- CREATE TABLE IF NOT EXISTS kmer (
--     kmer_id             serial_id UNIQUE          PRIMARY KEY,
--     kmer_sequence       VARCHAR(30) NOT NULL
-- );

-- CREATE TABLE IF NOT EXISTS sample_kmer (
--     sample_kmer_id      serial_id UNIQUE          PRIMARY KEY,
--     sample_id           varchar_id REFERENCES sample,
--     kmer_id             serial_id REFERENCES kmer,
--     sample_kmer_contig  serial_id REFERENCES contig,
--     sample_kmer_index   SMALLINT NOT NULL
-- );

CREATE TABLE IF NOT EXISTS kmer (
    seq         VARCHAR(30) NOT NULL,
    sample_id   VARCHAR(24) NOT NULL,
    contig      SMALLINT NOT NULL,
    index       SMALLINT NOT NULL,
    PRIMARY KEY(seq, sample_id, contig, index)
);

-- CREATE OR REPLACE FUNCTION kmer_upsert (
--     seq_         VARCHAR(30),
--     sample_id_   VARCHAR(24),
--     contig_      SMALLINT,
--     index_      SMALLINT
-- ) RETURNS VOID AS $$
-- BEGIN
--     INSERT INTO kmer 
--     ( seq , sample_id , contig , index  ) VALUES
--     ( seq_, sample_id_, contig_, index_ )
-- END;;
-- $$ LANGUAGE plpgsql;

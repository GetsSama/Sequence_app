import pandas as pd
import EntryModule as em


class Peptides:
    def __init__(self, path_to_sequence, path_to_ABL, transcrypt_name, drug_name):
        self.__sequence = em.Sequence_entity(path_to_sequence, transcrypt_name)
        self.__abl_data = em.ABL_entity(path_to_ABL, drug_name)
        self.__pept_tables = list()
        self.__is_tables_created = False

        for i in range(15):
            self.__pept_tables.append(list())

        self.__sequence.create_replaced_dict(self.__abl_data.replaced_positions, self.__abl_data.new_letters)

    def create_peptide_tables(self):
        if not self.__is_tables_created:
            drug = self.__abl_data.drug_name
            transcrypt = self.__sequence.transcrypt_name

            for pos in self.__abl_data.replaced_positions:
                left_part = ""
                right_part = ""
                pos = int(pos)
                mutation = self.__abl_data.get_mutation_by_position(str(pos))
                sequence_one_str = self.__sequence.get_sequence_by_position(str(pos))
                peptide = sequence_one_str[pos - 1]
                seq_len = len(sequence_one_str)

                if pos <= 15 or (pos + 15) > seq_len:
                    continue

                #counter += 1

                for i in range(1, 16):
                    left_part = sequence_one_str[pos - 1 - i]
                    right_part = sequence_one_str[pos - 1 + i]
                    peptide = left_part + peptide + right_part
                    self.__pept_tables[i - 1].append(list([peptide, drug, mutation, transcrypt]))

                self.__is_tables_created = True

    def peptides_to_csv(self):

        if self.__is_tables_created:
            df_list = list()

            for table in self.__pept_tables:
                df = pd.DataFrame(table, columns=['PEPTIDE', 'DRUG_NAME', 'AA_MUTATION', 'TRANSCRYPT'])
                df.index.name = 'NUMBERS'
                df_list.append(df)

            for i in range(len(df_list)):
                name = "pept_" + str(i + 1) + ".csv"
                df_list[i].to_csv(name)
                print("---> Create new file: " + str(name))

        else:
            print("Empty peptide tables, please create peptide tables and try again!")

    @property
    def peptide_tables(self):
        return self.__pept_tables
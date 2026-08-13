# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import tempfile
import unittest

import antismash
from antismash.common.secmet.locations import FeatureLocation
from antismash.common.secmet.test.helpers import DummyCDS, DummyRecord
from antismash.common.test.helpers import run_and_regenerate_results_for_module
from antismash.config import build_config, destroy_config, update_config
from antismash.detection import subclusters


DATA = {
    "AJAP_31990": {
        "location": [17866, 18657],
        "translation": (
            "MNGVRYEKKDHVAYVTMDRPEALNAMDRRMHAELAGIWDDVEADDDVRAVVLTGAGDRAFSVGQD"
            "LKERARLTDEGVEASTFGSSGQPGHPRLTDRFTLSKPVVARVHGYALGGGFELVLACDIVVASEE"
            "AVFGLPEVRLGLIAGAGGVFRLPRQLPQKVAMGYLLTGRRMDAATALRHGLVNEVVPFPELDRCV"
            "AEWTDDLVRAAPLSVRAIKEAALRSLDLPLEEAFKTSYPWEERRRTSADASEGARAFAEKRDPIW"
            "TGR"
        ),
    },
    "AJAP_31995": {
        "location": [18654, 19868],
        "translation": (
            "MTALSPGLDLRALAGEAHRVDDEVRTLRAAFVEAHAEELYAELTDGRTRYLRIDELVRAAALAFP"
            "GLVPSEAQMAAERARPQAEKEGREIDQGIFLRGILRAPKAGPHLLDAMLRPTARARRLLPEFLET"
            "GVVQMEAVRLERRDGVAHLTLCRDDCLNAEDAQQVDDMETAVDLVLLDPSVRVGLLRGGVMSHPR"
            "YLGRRVFCAGINLKKLSSGDIPLVDFLLRRETGYIQKIFRGLLTDDSWHSRFTGKPWMAAVDSFA"
            "IGGGTQLLLVFDHVLAASDAYLSLPAAKEGIIPGVSNFRLSRIAGPRVARQVILGGRKLRADEPD"
            "ARSIVDEVVPPEEMDAAIDGALARLDGEAVAANRRMVNLSEEPPEEFRRYIAEFALQQALRIYGA"
            "DVIGKVDGFAVGSR"
        )
    },
    "AJAP_32000": {
        "location": [18865, 20527],
        "translation": (
            "MVTRNEIKGELVLRFDGSRPLSAAAVEEIGAFCDRAEDQREPGPVTIHVTGAPPADWAKGLAVGL"
            "VSKWERVVRRFERLGRLTAAVASGECAGTALDLLLAADIRIAEPGTTLRLASAGGGTWPGMTVYR"
            "LTKQAGAAGIRRAVLLGAPIGTERALALDLIDEVSDDPAKTLAELDVVVDGAETAIRRQLIFEAG"
            "STTFEEALGSHLAAADRALRREAKS"
        )
    },
    "AJAP_32005": {
        "location": [20521, 21621],
        "translation": (
            "MTAIIEPAEDLSVLTGLTEITRFAGVGTAVSESSYSQTELLEILDIEDPKIRSVFLNSAIDRRFL"
            "TLPPEGPGGARATEPQGDLLDKHKRIAVDMGCRALEACLKSAGATLSDLRHLCCVTSTGFLTPGL"
            "SALIIREMGIDPHCSRSDIVGMGCNAGLNALNVVAGWSAAHPGELGVVLCSEACSAAYALDGTMR"
            "TAVVNSLFGDGSAALALVSGDGRVPGPRVLKFASYIITDAVDAMRYDWDRDQDRFSFFLDPQIPY"
            "VVGAHAEIVIDRLLSGTGLRRSDIGHWLVHSGGKKVIDAVVVNLGLSRHDVRHTTGVLRDYGNLS"
            "SGSFLFSYERLADEDVARPGDYGVLMTMGPGSTIEMALIQW"
        )        
    },
    "AJAP_32035": {
        "location": [27846, 28916],
        "translation": (
            "MTHLCLDDLERAARTALPGEIWDFLAGGSGAEASLAANRTALDRIFVIPRMLRELTGGTTEAEVL"
            "GRRASLPVAVAPVAYQRLFHPEGELAAARAARDAGVPYTICTLSSVPLEEIAVVGGRPWFQLYWL"
            "RDEKRSLELVRRAEDAGCEAIVFTVDVPWMGRRLRDMRNGFALPESVTAANFDAGAAAHRRTEGL"
            "SAVADHTAREFAPATWESVEAVRAHTNLPVVLKGILAVEDAVRAVDAGATGIVVSNHGGRQLDGA"
            "VPGIAMLEEIADAVSGGCEVLLDGGIRTGGDVLKALALGASGVLIGRPFMWGLAADGQAGARQVL"
            "DLLAVELRNALGLAGCDSVSAARRLSTRYER"
        )
    },
    "AJAP_32040": {
        "location": [28913, 29965],
        "translation" :(
            "MSSHENTQNFEIDYVEMYVANLEVAASGWLDKYDFSVTATDRSADHRSVTLRHSAIALVLTEPLS"
            "DRHPGATYLQTHGEGVADIALRTTDVAAAFEAAVKAGAKPLREPEKGVDSVLTATVSGFGDVVHT"
            "LIQSDVAEEAPRGKGGVELGVIDHFAVCLNAGDLGPTVAFYERALGFKQIFEEHIVVGAQAMNST"
            "VVQSTSGAVTLTLIEPDKTADPGQIDDFIKEHHGSGVQHIAFTSPDAVRAVKELSARGVEFLKTP"
            "DAYYDLLGERIELETHSLDDLRETKLLADEDHGGQLFQIFTASTHPRKTIFFEIIERQGAGTFGS"
            "SNIKALYEAVELERTGQSKLGPARR"
        )
    },
    "AJAP_32060": {
        "location": [33933, 35294],
        "translation" : (
            "MDNSRKRGVHRLFLTLRRTVEILVFMDSNGLSTHLNVETLHGSLTDPAISSMNLLNELIDEYPVA"
            "ISMAAGRPYEEFFDIRLIHTYLDAYCDHLRRDRKMDEAVVTRTLFQYGTTKGVIADLIAKNLAED"
            "ENIDAAPESVVVTVGAQEAMFLVLRTLRAGERDVLLAPAPTYVGLTGAALLTDTPVWPVRSNENG"
            "IDPDDLVLQLKRADEQGKRVRACYVTPNFANPTGTSMDLPSRHRLLDVAESNGILLLEDNAYGLF"
            "GSERLPSLKSLDRSGSVVYLGSFAKTGMPGARVGFAVADQRMADGGLFADQLSKLKGMLTVNTSP"
            "IAQAVIAGKLLLNDFSLTKANAREIAIYQRNLRLTLDALERGLGSCEGVSWNTPTGGFFVTVTVP"
            "FVVDDELLETAAREHGVLFTPMHHFYGGKGGFHQLRLSISLLTPELIEEGVARLAALIKPRLP"
        )
    },
    "AJAP_32155": {
        "location": [79880, 80710],
        "translation": (
            "MTIEKALVVGTGLIGTSAALSLREKGVAVHLSDIDAQAVRVARELGAGREWTGEEVDLAVIAVPP"
            "QLVGERLAELQKRGAARAYTDVASVKVDPIADAERLGCDMTSYVPGHPLAGRERSGPAAARAGLF"
            "AGRPWALCPGPGTGAEALRLTRGLVALCGAEAVTVGAAEHDSAVALVSHAPHVAASAVAASLASG"
            "DDVALSLAGQGLRDVTRIAAGNPLLWRRILSGNAVPVAAVLDRIAADLAAAATALRAGDLDDLTE"
            "LLRRGVDGHGRIPAIG"
        )
    }
}


class TestSubclusters(unittest.TestCase):
    def setUp(self):
        options = build_config(self.get_args(), isolated=True, modules=antismash.get_all_modules())
        self.options = update_config(options)
        subclusters.prepare_data()

    def tearDown(self):
        destroy_config()

    def get_args(self):
        return ["--minimal", "--subclusters", "--subclusters-region-mode", "create",
                "--enable-html"]

    def test_full_pathway(self):
        features = []
        for name, data in DATA.items():
            features.append(DummyCDS(locus_tag=name, start=data["location"][0],
                                     end=data["location"][1], translation=data["translation"]))
        record = DummyRecord(seq="A" * features[-1].end, features=features)
        with tempfile.NamedTemporaryFile(suffix=".gbk") as temp:
            record.to_genbank(temp.name)
            results = run_and_regenerate_results_for_module(temp.name, subclusters, self.options)

        assert results

        predicted_subregions = results.get_predicted_subregions()
        assert len(predicted_subregions) == 1
        assert predicted_subregions[0].label == "subclusters"
        assert predicted_subregions[0].location == FeatureLocation(17866, 80710, 1)
        assert predicted_subregions[0].tool == "subclusters"

        assert len(results.rule_results.protoclusters) == 2
        assert [(proto.product, proto.location) for proto in results.rule_results.protoclusters] == [
            ("SCG0041", FeatureLocation(17866, 35294, 1)),
            ("SCG0042", FeatureLocation(27846, 80710, 1)),
        ]

        assert len(results.predictions) == 2
        dhpg, hpg = results.predictions

        assert dhpg.rule.name == "SCG0041"
        assert dhpg.compound.name == "3,5-Dihydroxyphenylglycine (Dhpg)"
        assert dhpg.location == FeatureLocation(17866, 35294, 1)
        assert len(dhpg.domain_hits) == 6
        assert sorted((hit.cds_name, hit.domain_name) for hit in dhpg.domain_hits) == [
            ("AJAP_31990", "ECH_1"),
            ("AJAP_31995", "ECH_1"),
            ("AJAP_32000", "ECH_1"),
            ("AJAP_32005", "Chal_sti_synt_C"),
            ("AJAP_32005", "Chal_sti_synt_N"),
            ("AJAP_32060", "Aminotran_1_2"),
        ]
        aminotransferase = dhpg.domain_hits_by_cds["AJAP_32060"][0]
        assert aminotransferase.domain_accession == "PF00155.24"
        assert aminotransferase.domain_description == "Aminotransferase class I and II"

        assert hpg.rule.name == "SCG0042"
        assert hpg.compound.name == "4-Hydroxyphenylglycine (Hpg)"
        assert hpg.location == FeatureLocation(27846, 80710, 1)
        assert len(hpg.domain_hits) == 6
        assert sorted((hit.cds_name, hit.domain_name) for hit in hpg.domain_hits) == [
            ("AJAP_32035", "FMN_dh"),
            ("AJAP_32040", "Glyoxalase"),
            ("AJAP_32040", "Glyoxalase_4"),
            ("AJAP_32060", "Aminotran_1_2"),
            ("AJAP_32155", "PDH_C"),
            ("AJAP_32155", "PDH_N"),
        ]


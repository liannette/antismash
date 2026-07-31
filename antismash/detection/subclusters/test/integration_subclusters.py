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
        return ["--minimal", "--subclusters", "--subclusters", "any", "--enable-html"]

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
        assert len(results.rule_results.protoclusters) == 2

        assert [(proto.product, proto.location) for proto in results.rule_results.protoclusters] == [
            ("SCG0041", FeatureLocation(17866, 35294, 1)),
            ("SCG0042", FeatureLocation(27846, 80710, 1)),
        ]

        # assert len(results.get_predicted_subregions()) == 2
        # assert [(sub.label, sub.location, sub.tool) for sub in results.get_predicted_subregions()] == [
        #     ("SCG0041", FeatureLocation(17866, 35294, 1), "subclusters"),
        #     ("SCG0042", FeatureLocation(27846, 80710, 1), "subclusters"),
        # ]

        assert len(results.predictions) == 2
        dhpg, hpg = results.predictions

        assert dhpg.rule.name == "SCG0041"
        assert dhpg.location == FeatureLocation(17866, 35294, 1)
        # AJAP_32035 and AJAP_32040 fall within the cutoff of this rule,
        # but define the other subcluster, so they must not be included here
        assert len(dhpg.cds_results) == 5
        expected = {
            "AJAP_31990": ["ECH_1"],
            "AJAP_31995": ["ECH_1"],
            "AJAP_32000": ["ECH_1"],
            "AJAP_32005": ["Chal_sti_synt_C", "Chal_sti_synt_N"],
            "AJAP_32060": ["Aminotran_1_2"],
        }
        cds_mapping = {result.cds.get_name(): sorted(result.definition_domains[dhpg.rule.name])
                       for result in dhpg.cds_results}
        assert cds_mapping == expected

        aminotransferase = dhpg.domain_hits_by_cds["AJAP_32060"][0]
        assert aminotransferase.domain_accession == "PF00155.24"
        assert aminotransferase.domain_description == "Aminotransferase class I and II"

        assert hpg.rule.name == "SCG0042"
        assert hpg.core_location == FeatureLocation(27846, 80710, 1)
        assert len(hpg.cds_results) == 4
        assert {result.cds.get_name(): sorted(result.definition_domains["SCG0042"])
                for result in hpg.cds_results} == {
            "AJAP_32035": ["FMN_dh"],
            "AJAP_32040": ["Glyoxalase", "Glyoxalase_4"],
            "AJAP_32060": ["Aminotran_1_2"],
            "AJAP_32155": ["PDH_C", "PDH_N"],
        }

